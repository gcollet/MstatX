/* Copyright (c) 2010 Guillaume Collet
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE. 
 */

#ifndef __FACTORY_H_INCLUDED__
#define __FACTORY_H_INCLUDED__

#include <stdexcept>
#include <map>

/************************************
 * Factory
 ************************************/
template<class BaseClass>
class Factory
{
  public:
    typedef BaseClass* (*baseclass_returning_fun_type)();                    //a pointer to function returning pointer to BaseClass
    typedef std::map<std::string, baseclass_returning_fun_type> map_type;
    typedef typename map_type::iterator iterator;

  private:
    Factory(){};
    static
    map_type& Get()
    {
      static map_type map;
      return map;
    }

  public:
    static
    BaseClass*
    CreateByName(const std::string&name)
    {
      iterator it = Get().find(name);                       // look for the object
      if (it == Get().end())                                // object with name "name" is not found
        throw std::runtime_error(name + " not found\n");
      return it->second();									//Call the object Constructor... and return the object!
    }
  

    // This construct an object of type F and returns a pointer to it...
    template<class F>
    static
    BaseClass* Constructor()
    {
      return new F;
    }
    
    template<class F>
    static void Add(const std::string &name)
    {
      Get()[name] = &Constructor<F>; 
    }
};

#endif
