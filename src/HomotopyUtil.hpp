#include <vector>
#include <string>

#include <iostream>
#include <fstream>
#include <sstream>

#include <cstdio>

namespace Homotopy {

    namespace Util {
        namespace Matrix {

            std::vector< std::vector< double > > readCSV( const std::string name, const size_t rows ){
                
                std::vector< std::vector< double > > data;
                
                std::vector< std::string > entries;
                std::ifstream myfile (name.c_str());
                if (myfile.is_open())
                {
                    int l = 0;
                    std::string line;
                    while ( getline(myfile, line ) )
                    {
                        std::stringstream iss;
                        iss << line;
                        std::string token;
                        while ( std::getline (iss,token,',') ){
                            entries.push_back(token);
                        }
                    }
                    myfile.close();

                    //Infer the number of columns
                    size_t cols = entries.size() / rows;

                    for(size_t i=0; i<rows; i++){
                        std::vector< double > row;
                        for(size_t j=0; j<cols; j++){
                            double entry;
                            sscanf( entries[( i*cols + j )].c_str(), "%lf", &entry );
                            row.push_back( entry );
                        }
                        data.push_back( row );
                    }
                }

                else {
                    std::cout << "Unable to locate: " << name << std::endl;
                }

                return data;
            }

            void SplitAtCol( const std::vector< std::vector< double > >& data, const size_t loc, std::vector< std::vector< double > >& left, std::vector< std::vector< double > >& right ){
                std::vector< std::vector< double > >::const_iterator i;
                for( i=data.begin(); i!=data.end(); i++){
                    std::vector< double > row_left( i->begin(), i->begin() + loc );
                    left.push_back( row_left );

                    std::vector< double > row_right( i->begin() + loc, i->end() );
                    right.push_back( row_right );
                }
            }

            std::vector< std::vector< double > > GetEye( const size_t dim ){
                std::vector< std::vector< double > > ret;
                /*std::vector< double > row( dim, 0.0 );
                for( size_t i=0; i<dim; i++ ){
                    ret.push_back( row );
                    ret.back()[i] = 1.0;
                }*/
                std::vector< double > o(3, 0.0);
                o[0] = 1.0; o[1] = 1.0 ; o[2] = 0.0;
                ret.push_back( o );
                o[0] = 1.0; o[1] = 0.0 ; o[2] = 1.0;
                ret.push_back( o );
                o[0] = 0.0; o[1] = 1.0 ; o[2] = 1.0;
                ret.push_back( o );
                return ret;
            }
        }
    }
}