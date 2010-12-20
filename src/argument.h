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
		virtual void setValue(string val) = 0;
		
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
