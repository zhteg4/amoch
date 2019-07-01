%{
#include <cstdio>
#include "parse.hh"
%}

%option noyywrap yylineno nounput

%%
#.*\n                 // ignore comment to end of line
\[                    {return TOK_LBRACKET;}
\]                    {return TOK_RBRACKET;}
,                     {return TOK_COMMA;}
\{                    {return TOK_LBRACE;}
\}                    {return TOK_RBRACE;}
\'[^'\n]*\'           {snprintf(yylval.s, 128, "%s", yytext); return TOK_QSTR;}
\"[^"\n]*\"           {snprintf(yylval.s, 128, "%s", yytext); return TOK_QSTR;}
[^ \t\n\{\}\[\],]+    {snprintf(yylval.s, 128, "%s", yytext); return TOK_STR;}
[ \t\n]+              // ignore whitespace
%%
