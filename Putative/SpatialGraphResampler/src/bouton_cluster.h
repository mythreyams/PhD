#include "bouton_params.h"
#include "stdlib.h"

class bouton_cluster{
    
public:
    bouton_cluster(){};
    bouton_cluster(int, std::list<double *>*,std::list<double *>* );
    std::list<double *>* get_auto_boutons(void)
    {
        return auto_boutons;
        
    };
     
private:

    // data members
    enum{
        Init,
        Compare,
        Push,
        End
    }state;
    std::list<double *> *auto_boutons;
    std::list<double *> rawBoutonList;
    std::list< std::list<double*> *> groupedBoutonLists;
    
    // private function
    void getInitialCentroids( int, double *centroid1, double *centroid2);
    void assignClusters( int, double *centroid1, double *centroid2, bool *done_with_clustering);
    void computeNewCentriods( int, double *centroid1, double *centroid2);
    void createRowBoutonsList( int NumOfLandmarks);
    void groupBoutonPoints( void);
    void findBoutonandOverlapCentroids( int NumOfLandmarks, std::list<double *>* bouton_centroids, std::list<double *>* overlap_centroids);
    
};