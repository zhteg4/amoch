%{
// ============================================================================
//
// parse_input() is called BEFORE any log file is open, so write error messages
// to stderr 
//

#include <cstdio>
#include <cstring>
#include <cerrno>
#include "param.h"

using AMOCH::Param;

static AMOCH::ParamParser parser;
static const char *yyfile = NULL;
static std::vector<std::string> pval(0, "");

// Defined in lex.cc
extern int   yylineno;
extern char *yytext;
extern FILE *yyin;
extern int   yylex();
extern void  yyrestart(FILE *);
extern int   yylex_destroy(void);

// ============================================================================
// Parse error message
void
yyerror(const char *msg)
{
   fprintf(stderr, "Parse error at \"%s\" in input file %s, line %d: %s\n", 
           yytext, yyfile, yylineno, msg);
}
           
// ============================================================================
// Read input file into Param tree with root proot; return true on success, 
// false on error
// ============================================================================
int yyparse();

bool 
parse_input(Param *proot,
            const char *path)
{
   bool ret = true;
   const char *proot_name = proot->key().c_str();

   if (!(yyin = fopen(path, "r"))) {
      fprintf(stderr, "Unable to open input file %s: %s\n", path, 
              strerror(errno));
      return false;
   }
   yyfile = path;
   yyrestart(yyin);
   parser.push(proot);
   if (0 != yyparse())
      ret = false;
   fclose(yyin);
   yylex_destroy();

   const char *name = parser.stack_key();

   if (name != proot_name) {
      fprintf(stderr, "Parse error: value \"%s\" left on stack\n", name);
      ret = false;
   }
   (void)parser.pop();  // remove proot -- assumes proot is cleaned up elsewhere
   return ret;
}
%}

%union {
   char s[128];
}

%token TOK_LBRACKET TOK_RBRACKET TOK_COMMA TOK_LBRACE TOK_RBRACE 
%token <s> TOK_STR TOK_QSTR

%%
input
   : keys
   | input keys
   ;

keys
   : keyval
   | keyvalbr
   | endbr
   ;

items
   : TOK_STR                  {pval.push_back($1);}
   | items TOK_COMMA TOK_STR  {pval.push_back($3);}
   ;

list: TOK_LBRACKET items TOK_RBRACKET
   ;

keyval
   : TOK_STR TOK_STR
   {
      Param *p = new Param($1);

      parser.add_value(p, $2);
      parser.add_sub(p);
   }
   | TOK_STR TOK_QSTR
   {
      Param *p = new Param($1);

      std::string qstr($2);
      qstr.erase(qstr.begin()+qstr.size()-1);
      qstr.erase(qstr.begin());
      parser.add_value(p, qstr);
      parser.add_sub(p);
   }
   | TOK_STR list
   {
      Param *p = new Param($1);

      parser.set_values(p, pval);
      pval.clear();
      parser.add_sub(p);
   }
   ;

keyvalbr: TOK_STR TOK_STR TOK_LBRACE
   {
      Param *p = new Param($1);
      
      parser.add_value(p, $2);
      parser.push(p);
   }
   ;

endbr: TOK_RBRACE
   {
      if (parser.stack_key() == NULL) {
         yyerror("mismatched braces");
         YYABORT;
      }
      else
         parser.add_sub(parser.pop());
   }
   ;
%%
