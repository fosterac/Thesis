#include <vector>
#include <string>

#include <iostream>
#include <fstream>

namespace Homotopy {

    namespace Util {
        namespace Matrix {

            std::vector< std::vector< double > > readCSV( const std::string name, const size_t rows, const size_t cols ){
                std::vector< std::string > entries;
                std::ifstream myfile (name.c_str());
                if (myfile.is_open())
                {
                    while ( myfile.good() )
                    {
                        std::string line;
                        std::getline (myfile,line,',');
                        entries.push_back(line);
                    }
                myfile.close();
                }

                std::vector< std::vector< double > > data;
                for(size_t i=0; i<rows; i++){
                    std::vector< double > row;
                    for(size_t j=0; j<cols; j++){
                        row.push_back( entries.at( i*cols + j ) );
                    }
                    data.push_back( row );
                }

                return data;
            }

            void SplitAtCol( const std::vector< std::vector< double > >& data, const size_t loc, std::vector< std::vector< double > >& left, std::vector< std::vector< double > >& right ){
                std::vector< std::vector< double > >::iterator i;
                for( i=data.begin(); i!=data.end(); i++){
                    std::vector< double > row_left( i->begin(), i->begin() + loc );
                    left.push_back( row_left );

                    std::vector< double > row_right( i->begin() + loc, i->end() );
                    right.push_back( row_right );
                }
            }

            std::vector< std::vector< double > > GetEye( const size_t dim ){
                std::vector< std::vector< double > > ret;
                std::vector< double > row( dim, 0.0 );
                for( size_t i=0; i<dim; i++ ){
                    ret.push_back( row );
                    ret.back()[i] = 1.0;
                }
                return ret;
            }
        }
    }
}