// ============================================================================
// param.h -- AMOCH::Param, AMOCH::ParamParser classes
// ----------------------------------------------------------------------------
// Copyright
// ----------------------------------------------------------------------------
// LICENSE
// ============================================================================

#ifndef AMOCH_PARAM_H
#define AMOCH_PARAM_H

#include <string>
#include <vector>

namespace AMOCH {

class ParamParser;

// ----------------------------------------------------------------------------
//
// Param -- input parameter key - value 
//
class Param {

   friend class ParamParser;

   private:

   std::string _key;
   std::vector<std::string> _value;
   Param *_sub;   // down a level: level1 { level2 ... }
   Param *_next;  // next Param at current level

   // Disable copy constructor and assignment
   Param(const Param&);
   Param& operator=(const Param&);

   // Return a pointer to the first Param in _next with matching key, or NULL
   // if no match
   const Param *_find(const std::string& key) const;

   public:

   explicit Param(const std::string& key);

   ~Param();

   // Return the key
   const std::string& key() const {return _key;}

   // Return the number of values
   int nvalues() const {return _value.size();}

   // Return a pointer to the nth value
   const std::string& value(int n=0) const {return _value[n];}

   // Return the first Param in _sub with matching key, or NULL if no match 
   const Param *find(const std::string& key) const
   {
      return ((_sub)) ? _sub->_find(key) : NULL;
   }

   // Return the next Param in _next with matching key, or NULL if no matches
   const Param *find_next(const std::string& key) const
   {
      return ((_next)) ? _next->_find(key) : NULL;
   }

};  // Param

// ----------------------------------------------------------------------------
//
// Used by the parser (parse.y) to manipulate Param objects while parsing input.
// Not used outside the parser.
//
class ParamParser {

   private:

   Param *_pstack;

   // Disable copy constructor and assignment
   ParamParser(const ParamParser&);
   ParamParser& operator=(const ParamParser&);

   public:

   ParamParser() : _pstack(NULL) {}

   ~ParamParser();
   
   // Return key of Param on stack; return NULL if stack is empty
   const char *stack_key() const
   {
      return ((_pstack)) ? _pstack->key().c_str() : NULL;
   }

   // Add a new value string to p
   void add_value(Param *p,
                  const std::string& val)  
   {
      p->_value.push_back(val);
   }

   // Add multiple value strings to p
   void set_values(Param *p,
                   const std::vector<std::string>& v)  
   {
      p->_value = v;
   }

   // Add p to the end of _pstack->sub
   void add_sub(Param *p);
      
   // Push p onto stack
   void push(Param *p);

   // Return Param popped from stack; return NULL if stack empty
   Param *pop();

};  // ParamParser

}   // AMOCH

#endif  // AMOCH_PARAM_H
