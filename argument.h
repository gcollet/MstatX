/*
 *  argument.h
 *  ARGU
 *
 *  Created by Guillaume Collet on 21/08/10.
 *
 */

#ifndef _ARGUMENT_H_
#define _ARGUMENT_H_

#include<string>
#include<exception>
#include <sstream>

using namespace std;

namespace ARGU {

	/**
	 * This class represents any argument
	 * (integer, float, string). Thus, there is 
	 * only one virtual function : setValue
	 */
	class Arg
	{
	protected:
		string _small_flag;		/**< The small version of the argument flag */
		string _long_flag;		/**< The long version of the argument flag */
		string _description;	/**< The description of the argument (for usage printing)*/
		string _type;					/**< The type of the argument (integer, float, string) */
		bool _needValue;			/**< This boolean value is true when a value is needed */
		
		Arg(const string & small_flag,
				const string & long_flag,
				const string & description,
				const string & type):
		_small_flag(small_flag),
		_long_flag(long_flag),
		_description(description),
		_type(type)
		{};
		
	public:
		/**
		 * The destructor is virtual
		 */
		virtual ~Arg(){};
		
		/**
		 * Returns the argument small flag.
		 */
		const string& getSmallFlag() const {return _small_flag;};
		
		/**
		 * Returns the argument long flag.
		 */
		const string& getLongFlag() const {return _long_flag;};
		
		/**
		 * Returns the argument description.
		 */
		string getDescription() const {return _description;};
		
		/**
		 * Returns the argument type.
		 */
		string getType() const {return _type;};
		
		/**
		 * Indicates whether the argument needs a value or not.
		 */
		bool needValue() const {return _needValue;};
		
		/**
		 * Sets the value of the argument.
		 */
		virtual void setValue(string val){};
		
	};
	
	/**
	 * ValueArg is an argument with a value of type T
	 */
	template <typename T>
	class ValueArg : public Arg {
	private:
		T _value;
	
	public:
		ValueArg(const string & small_flag,
						 const string & long_flag,
						 const string & description,
						 const string & type,
						 const T value):
		Arg(small_flag,long_flag,description,type),
		_value(value)
		{
			_needValue = true;
		};
		
		~ValueArg(){};
		
		T getValue() const {return _value;};
		
		void setValue(string val) {
		  istringstream is(val);
			while (is.good()) {
				is >> _value;
			}
		};
	};
	
	class SwitchArg : public Arg {
	private:
		bool _value;
	public:
		SwitchArg(const string & small_flag,
						  const string & long_flag,
						  const string & description,
							const string & type):
		Arg(small_flag,long_flag,description,type),
		_value(false)
		{
			_needValue = false;
		};
		
		~SwitchArg(){};
		
		bool getValue() const {return _value;};
		void setValue(string val) {_value = true;};
	};
}

#endif