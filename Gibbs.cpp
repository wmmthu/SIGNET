#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector Gibbs(SEXP pNeighbor_statistics, SEXP pNeighbor_index, SEXP pNeighbor_weight, SEXP pNeighbor_sumW, SEXP ppara, SEXP pp_value, SEXP pZ)
{
	//IntegerVector Gibbs(NumericMatrix Neighbor_statistics, List& Neighbor_index, List& Neighbor_weight, NumericVector& para, NumericVector& p_value, IntegerVector& Z)

	NumericMatrix Neighbor_statistics(pNeighbor_statistics);
	List Neighbor_index(pNeighbor_index);
	List Neighbor_weight(pNeighbor_weight);
	NumericMatrix Neighbor_sumW(pNeighbor_sumW);
	NumericVector para(ppara);
	NumericVector p_value(pp_value);
	IntegerVector Z(pZ);

	int num = p_value.size();
	int nrow = Neighbor_statistics.nrow();
	int ncol = Neighbor_statistics.ncol();

	if(num != nrow){
		printf("Error : number of genes != nrow of neighbor statistics matrix.\n");
	}

	for(int g=0;g<num;g++){

		double tmp = 0;
		for(int j=0;j<ncol;j++){ tmp+= Neighbor_statistics(g,j)*para(j+3);}
		int before = Z(g);

		//tmp = exp(tmp+para(2)) * para(1)*pow(p_value(g),para(1)-1)/(para(0)*pow(p_value(g),para(0)-1));
		tmp = tmp + para(2) + log(para(1)/para(0)) + (para(1)-para(0))*log(p_value(g));
		//Z(g) = rbinom(1,1, tmp/(1+tmp))[0];
		Z(g) = rbinom(1,1, 1/(1+exp(-tmp) ))[0];
		
		int after = Z(g);
		
		if(before == after) continue;
		for(int k=0;k<ncol;k++){
			IntegerVector neighbors(as<List>(Neighbor_index[k])(g));
			NumericVector tmpW(as<List>(Neighbor_weight[k])(g));

			for(int t=0;t<neighbors.size();t++){
				
				int i = neighbors(t) - 1;
				double w = tmpW(t);
				//Neighbor_statistics(i,k) = Neighbor_statistics(i,k)  + (after-before) * 2 * w / Neighbor_sumW(k,i);
				Neighbor_statistics(i,k) = Neighbor_statistics(i,k)  + (after-before) * 2 * w;
				//Neighbor_statistics(i,k) = Neighbor_statistics(i,k)  + (after-before) * w;
				
			}
		}
		
		
	}
	return(Z);
}