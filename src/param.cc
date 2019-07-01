// ============================================================================
// param.cc -- AMOCH::Param methods
// ----------------------------------------------------------------------------
// Copyright
// ----------------------------------------------------------------------------
// LICENSE
// ============================================================================

#include "param.h"

using std::string;
using AMOCH::Param;
using AMOCH::ParamParser;

// ============================================================================
// ParamParser destructor
// ============================================================================
ParamParser::~ParamParser()
{
   Param *p = _pstack;

   while ((_pstack)) {
      p = _pstack->_next;
      delete _pstack;
      _pstack = p;
   }
}
   
// ============================================================================
// Add p to the end of _pstack->sub
// ============================================================================
void 
ParamParser::add_sub(Param *p)
{
   Param **pp = &_pstack->_sub;
      
   while ((*pp)) 
      pp = &(*pp)->_next;
   *pp = p;
}

// ============================================================================
// Push p onto stack
// ============================================================================
void 
ParamParser::push(Param *p)
{
   p->_next = _pstack;
   _pstack = p;
}

// ============================================================================
// Return Param popped from stack; return NULL if stack empty
// ============================================================================
Param *
ParamParser::pop()
{
   Param *p = _pstack;

   if ((p)) {
      _pstack = p->_next;
      p->_next = NULL;
   }
   return p;
}

// ============================================================================
// Param constructor
// ============================================================================
Param::Param(const string& key) :
   _key(key),
   _value(0, ""),
   _sub(NULL),
   _next(NULL)
{}

// ============================================================================
// Param destructor
// ============================================================================
Param::~Param()
{
   if ((_next)) {
      delete _next;
      _next = NULL;
   }
   if ((_sub)) {
      delete _sub;
      _sub = NULL;
   }
}

// ============================================================================
// Return a pointer to the first Param in _next with matching key, or NULL
// if no match
// ============================================================================
const Param *
Param::_find(const std::string& key) const 
{
   if (key == _key)
      return this;
   else if ((_next))
      return _next->_find(key);
   else 
      return NULL;
}

