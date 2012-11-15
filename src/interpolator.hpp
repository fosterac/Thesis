#ifndef interpolator_hpp
#define interpolator_hpp

#include "stdafx.h"
#include "interpolation.h"

#include <vector>
#include <algorithm>
#include <stdio.h>

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <string>

namespace Interpolation {
	class Interpolator { };

	typedef std::vector< std::vector< double > > Data;

	class RBF_AlgLib {
	private:
		int Points;
		int dim;
		
		std::vector< double > mins;
		std::vector< double > ranges;

		alglib::rbfmodel model;
		alglib::rbfreport report;

		template< typename T >
		class vertical_iterator : public T::const_iterator {
			int column;
		public:
			vertical_iterator( int c ) : column( c ) {}
			typename T::value_type::value_type operator*() const { return (this->T::const_iterator::operator*())[column]; }
			vertical_iterator& operator= (const typename T::const_iterator & rhs)  { this->T::const_iterator::operator=(rhs); return *this; }
		};

		static std::vector< double > FindMins(const Data &d ){
			//Make sure we have some data
			assert ( !d.empty() );
			assert ( !d.front().empty() );

			std::vector< double > result( d[0].size() );
			int i;
			for(i=0; i<result.size(); i++) { 
				vertical_iterator< Data > vi1(i);
				vi1 = d.begin();
				vertical_iterator< Data > vi2(i);
				vi2 = d.end();
				result[i] = *std::min_element( vi1, vi2 ); 
			}
			return result;
		}
		//Maybe we can eliminate dimensions that have no range...
		static std::vector< double > FindRanges(const Data &d ){
			std::vector< double > result( d[0].size() );
			int i;
			for(i=0; i<result.size(); i++) { 
				vertical_iterator< Data > vi1(i);
				vertical_iterator< Data > vi2(i);
				vi1 = d.begin();
				vi2 = d.end();
				result[i] = *std::max_element( vi1, vi2 ) - *std::min_element( vi1, vi2 ); 
				//assert( result[i] != 0 );
			}
			return result;
		}
		inline static double Rescale( double x, double min, double range) { return (range == 0) ? x : ((x-min)/range); }
		inline static double Unscale( double x, double min, double range) { return (range == 0) ? x : ((x*range)+min); }
		double* PrepareData( const Data &data, int N, int dim ) {
			double* d;
			try {
				d = new double[N*dim];
			}
			//catch( std::bad_alloc& ) {
			catch( std::exception& ) {
				printf("Dynamic allocation failed while staging RBF data for copy\n");
				d = NULL;
				return d;
			}

			int i;
			for(i=0; i<N; i++) {
				int j;
				for(j=0; j<dim; j++) { d[i*dim + j] = Rescale( data[i][j], this->mins[j], this->ranges[j] ); }
			}

			return d;
		}

	public:
		RBF_AlgLib( Data& data ) : mins( FindMins( data ) ), ranges( FindRanges( data ) ) {
			try
			{
				//Assume one dependent variable
				this->dim = data.front().size() - 1;
				this->Points = data.size();

				assert( this->dim == 2 || this->dim == 3 );

				double *d = PrepareData( data, this->Points, this->dim + 1);

				alglib::real_2d_array vals;
				vals.setcontent(this->Points, this->dim+1, d);
				delete [] d;

				alglib::rbfcreate(this->dim, 1, this->model);
				alglib::rbfsetpoints(this->model, vals);
				//alglib::rbfsetalgoqnn(this->model);
                rbfsetalgomultilayer(this->model, 1.0, 5, 1.0e-2);
				alglib::rbfbuildmodel(this->model, this->report);
				assert( int(this->report.terminationtype) == 1 );
			}
			catch(alglib::ap_error e)
			{
				printf("error msg: %s\n", e.msg.c_str());
			}			
		}

		double evaluate( const std::vector< double > &x){
			assert( x.size() == this->dim );
			if( this->dim == 2 ) 
				/*return Unscale( alglib::rbfcalc2(this->model, 
                Rescale(x[0], this->mins[0], this->ranges[0] ), 
				Rescale(x[1], this->mins[1], this->ranges[1] ) ), 
				this->mins[2], this->ranges[2] ) ;*/
                /*return alglib::rbfcalc2(this->model, 
				Rescale(x[0], this->mins[0], this->ranges[0] ), 
				Rescale(x[1], this->mins[1], this->ranges[1] ) );*/
                return alglib::rbfcalc2(this->model, 
				x[0], 
				x[1]);

			if( this->dim == 3 ) 
				return Unscale( alglib::rbfcalc3(this->model, 
				Rescale(x[0], this->mins[0], this->ranges[0] ), 
				Rescale(x[1], this->mins[1], this->ranges[1] ),
				Rescale(x[2], this->mins[2], this->ranges[2] ) ),
				this->mins[3], this->ranges[3] ) ;
		}
	};

	typedef RBF_AlgLib RBF;

    Data GetDataFromFile(const char* filename){
        Data ret;
        std::string line;
        std::ifstream myfile (filename);
        
        if (myfile.is_open())
        {
            while ( myfile.good() )
            {
                std::getline (myfile,line);

                std::vector< double > v;
                double x, y, z;

                sscanf( line.c_str(), "%lf %lf %lf", &x, &y, &z );

                v.push_back(x); v.push_back(y); v.push_back(z);
                ret.push_back(v);
            }
            myfile.close();
        }

        else throw std::runtime_error("Unable to open data file!");

        return ret;
    }
}

#endif