/* Boost.Flyweight example of use of key-value flyweights.
 *
 * Copyright 2006-2008 Joaquin M Lopez Munoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See http://www.boost.org/libs/flyweight for library home page.
 */

#include <boost/array.hpp>
#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using namespace boost::flyweights;

class ligand
{
public:
  ligand(const std::string& filename):filename(filename)
  {
    std::cout<<"loaded "<<filename<<" file"<<std::endl;
  }

  ligand(const ligand& x):filename(x.filename)
  {
    std::cout<<"ligand["<<filename<<"] copied"<<std::endl;
  }

  ~ligand()
  {
    std::cout<<"unloaded "<<filename<<std::endl;
  }

  const std::string& get_filename()const{return filename;}

  size_t get_value() const { return 5; }
  std::vector<size_t> atoms;

private:
  std::string filename;
};

struct ligand_filename_extractor
{
  const std::string& operator()(const ligand& x)const
  {
    return x.get_filename();
  }
};

typedef flyweight<key_value<std::string,ligand,ligand_filename_extractor> > ligand_flyweight;

int main()
{
  const char* ligand_filenames[]={"frag1.mol2", "frag2.mol2", "frag3.mol2", "frag4.mol2"};
  const int num_ligand_filenames=sizeof(ligand_filenames)/sizeof(ligand_filenames[0]);

  for(int i=0;i<50000;++i){
    ligand_flyweight lig(ligand_filenames[std::rand()%num_ligand_filenames]);
    size_t a = lig().get_value();
  }

  return 0;
}
