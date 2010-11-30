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
