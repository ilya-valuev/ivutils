#ifndef NAMED_FORMATER_H
#define NAMED_FORMATER_H

#include <string>
#include <vector>
#include <iostream>
#include <set>
//#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/function.hpp>
//#include <boost/variant.hpp>
#include <boost/bind.hpp>

/**
  Usage:
  \code
    // Simple case:
        std::cout
            <<  named_formater(" User: '$(user)' at $IP been locked.")
                              ("IP", "127.0.0.0")
                              ("user", "I'am")
            <<std::endl;
  \endcode
  \code
    // Use same set of key(s) at multiple templates
     named_formater::args argumens= named_formater::mk_arguments()
        ("user", "I'am")
        ("IP",   "127.0.0.0")
        ("more information",   "<<nothing to say>>")
     ;

     std::cout<<named_formater(" User: '$(user)' at $IP been locked.")(argumens)<<std::endl;
     std::cout<<named_formater(" Are '$(user)' same as 'i-am'?")(argumens)<<std::endl;
  \endcode
**/
class named_formater;


/**
  This function can be used at case of time-expansive viarble caluation.
  Usage:
  \code

    std::string A();
    std::string B();
    std::string encoder(const std::string& name)
    {
        if(name=="A") return A();
        if(name=="B") return B();
        // ...
        return "Wrng key...";
    }
    //...
    std::cout<<named_format<std::string>(" $A+$B ", '$', &encoder )<<std::endl;
  \endcode
**/
template< typename StringType, typename ValueByNameFunction>
StringType named_format(const    StringType&            template_string,
                        typename StringType::value_type sign,
                        const    ValueByNameFunction&   value_by_name);



class named_formater
{
    static std::string throw_on_no_argument(std::string name);
public:
    struct arg
    {
        std::string name, value;
        bool recursive;

        arg( const std::string& name)
            : name(name), value(),  recursive(false)
        {}

        arg( const std::string& name, const std::string& value, bool recursive=false)
            : name(name), value(value), recursive(recursive)
        {}

        template <typename T>
        arg( const std::string& name, const T& value , bool recursive=false)
            : name(name), value( boost::lexical_cast<std::string>(value) ), recursive(recursive)
        {}


        bool operator == (const arg& rhs) const { return name == rhs.name; }
        bool operator <  (const arg& rhs) const { return name <  rhs.name; }
    };
    struct args : std::set<arg>
    {
        args& add(  const arg& a)
        {
            erase(a.name);
            insert(a);
            return *this;
        }
        template <typename T>
        args& operator() (  const std::string& name, const T& value, bool recursive=false )
        {
            return add(arg(name, value, recursive));
        }
    };

    static args mk_arguments() { return args(); }

    named_formater(const std::string& _t, char sign='$' ) : m_template_sign(sign), m_template(_t){}

    typedef std::string f_on_no_argument(std::string name);
    const std::string& str(f_on_no_argument *on_no_argument=throw_on_no_argument) const;
    operator std::string () const { return str(); }



    named_formater operator () ( const arg& a ) const
    {
        named_formater ret(*this);
        ret.m_args.add( a );
        return ret;
    }
    named_formater operator () ( const args& a ) const
    {
        named_formater ret(*this);
        ret.m_args.insert( a.begin(), a.end() );
        return ret;
    }

    template <typename T>
    named_formater operator () ( const std::string& name, const T& value, bool recursive=false ) const
    {
        named_formater ret(*this);
        ret.m_args.add( arg(name, value, recursive) );
        return ret;
    }
private:
    mutable std::string m_ret;
    char        m_template_sign;
    std::string m_template;
    args        m_args;

    std::string value_by_name( const std::string & v, f_on_no_argument* on_no_argument ) const;
};

inline std::ostream& operator << ( std::ostream& out, const named_formater & fmt )
{
    return out<<fmt.str();
}

inline std::string named_formater::throw_on_no_argument(std::string name)
{
    throw std::runtime_error(
                named_formater("Argument '$(name)' not found.")
                              ("name", name)
                              .str().c_str()
                );
}

template< typename StringType, typename ValueByNameFunction>
StringType named_format(const    StringType&            template_string,
                        typename StringType::value_type template_sign,
                        const    ValueByNameFunction&   value_by_name)
{
    StringType ret;
    typedef typename StringType::size_type t_pos;
    t_pos i;
    for(i=0; i!=StringType::npos && i<template_string.size();  )
    {
        t_pos j=template_string.find(template_sign, i);

        if(j!=0 && j<template_string.size() && template_string[j-1]=='\\')
        {
            ret += template_string.substr(i, j-i-1);
            ret+="$";
            i=j+1;
            continue;
        }
        else if(j==StringType::npos || j==template_string.size())
        {
            ret += template_string.substr(i);
            break;
        }
        else
        {
            ret += template_string.substr(i, j-i);
        }

        StringType name;
        t_pos j_end;
        if(j+1==template_string.size())
        {
            name = "";
            j_end = j+1;
        }
        else
        {
            if(template_string[j+1]=='(')
            {
                j_end=template_string.find(')', j+1);
                if(j_end==std::string::npos)
                    j_end=template_string.size()+1;
                else
                    j_end=j_end+1;
                name = template_string.substr(j+2, j_end-j-3);
            }
            else
            {
                j_end=j+1;
                while( j_end<template_string.size() && ( isalnum( template_string[j_end] ) || template_string[j_end]=='_' ) )
                    ++j_end;
                name = template_string.substr(j+1, j_end-j-1);
            }
        }

        ret += value_by_name(name);
        i=j_end;
    }

    return ret;
}

inline const std::string& named_formater::str(f_on_no_argument* on_no_argument) const
{
    m_ret=named_format( m_template, m_template_sign, boost::bind( & named_formater::value_by_name ,this, _1, on_no_argument ) );
    return m_ret;
}

inline std::string named_formater::value_by_name(const std::string &name, named_formater::f_on_no_argument *on_no_argument) const
{
    args::const_iterator i= m_args.find(arg(name));
    if(i!=m_args.end())
    {
        if(i->name==name)
        {
            if(i->recursive)
                return named_format(i->value,
                                    m_template_sign,
                                    boost::bind( &named_formater::value_by_name, this, _1, on_no_argument ) );
            else
                return i->value;
        }
    }
    return on_no_argument(name);
}

#ifdef NAMED_FORMAT_TEST_CASESS
    bool test_named_formater_1()
    {
        return  named_formater("$test1+$(test2) \\$(test2)")
                ("test1", 1)
                ("test2", 2)
                .str() == "1+2 $(test2)";
    }
    bool test_named_formater_2()
    {
        return  named_formater("+$test1$(test2) \\$(test2)=\\$$test2")
                    ("test1", 1)
                    ("test2", 2)
                .str() == "+12 $(test2)=$2";
    }
    bool test_named_formater_3()
    {
        return  named_formater("$test1==2")
                ("test1", "$test2", true)
                ("test2", 2)
                .str() == "2==2";
    }
    bool test_named_formater_4()
    {
        return  named_formater("$A==2")
                ("A", 1)
                ("B", 2)
                ("A", 2)
                .str() == "2==2";
    }
#endif

#endif // NAMED_FORMATER_H
