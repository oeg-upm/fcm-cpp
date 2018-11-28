#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <float.h>
#include "fcm.h"


FCM::FCM(double m, double epsilon){
    m_epsilon = epsilon;
    m_m = m;
    m_membership = nullptr;
    m_data = nullptr;
    m_cluster_center = nullptr;
    m_num_clusters = 0;
    m_num_dimensions = 0;
}

FCM::~FCM(){
    if(m_data!=nullptr){
        delete m_data;
        m_data = nullptr;
    }

    if(m_membership!=nullptr){
        delete m_membership;
        m_membership = nullptr;
    }

    if(m_cluster_center!=nullptr){
        delete m_cluster_center;
        m_cluster_center = nullptr;
    }
}


double FCM::update_membership(){
    /*
     *
    */
    long k, i;
    double new_uik;
    double max_diff = 0.0, diff;

    if(m_data==nullptr || m_data->rows()==0){
        throw std::logic_error("ERROR: data should not be empty when updating the membership");
    }

    if(m_membership==nullptr || m_membership->rows() == 0 || m_membership->rows() != m_data->rows()){
        //cout << "init the membership";
        this->init_membership();
    }
    if(m_num_clusters==0){
        throw std::logic_error("ERROR: the number of clusters should be set");
    }

//    cout <<"mdata rows: "<< m_data->rows()<<endl;
//    cout << "mdata cols: "<<m_data->cols()<<endl;

    for (i = 0; i < m_num_clusters; i++) {
        for (k = 0; k < m_data->rows(); k++) {
            //cout << "point: " << k << " and cluster" << i <<endl;
            //cout << "\nwill ask for the new new_uik"<< endl;
            new_uik = this->compute_membership_point(i, k);
            //cout << new_uik << endl;
            diff = new_uik - (*m_membership)(k,i); // We need the membership inversed which is more natural for us
            if (diff > max_diff){
                max_diff = diff;
            }
            (*m_membership)(k,i) = new_uik;
        }
    }
    return max_diff;
}


void FCM::compute_centers(){
    long i, j, k;
    double numerator, denominator;
    MatrixXf t;
    t.resize(m_data->rows(), m_num_clusters);
    if(m_data == nullptr || m_data->rows() == 0){
        throw std::logic_error("ERROR: number of rows is zero");
        return;
    }
    for (i = 0; i < m_data->rows(); i++) { // compute (u^m) for each cluster for each point
        for (j = 0; j < m_num_clusters; j++) {
            t(i,j) = pow((*m_membership)(i,j), m_m);
        }
    }
    for (j = 0; j < m_num_clusters; j++) { // loop for each cluster
        for (k = 0; k < m_num_dimensions; k++) { // for each dimension
            numerator = 0.0;
            denominator = 0.0;
            for (i = 0; i < m_data->rows(); i++) {
                numerator += t(i,j) * (*m_data)(i,k);
                denominator += t(i,j);
            }
            (*m_cluster_center)(j,k) = numerator / denominator;
        }
    }
}

double FCM::get_dist(long i, long k){
  /*
   * distance which is denoted in the paper as d
   * k is the data point
   * i is the cluster center point
  */
  //cout<<"get_dist: point: "<<k<<" and cluster "<<i<<endl;
  long j;
  double sqsum = 0.0;
  if(m_num_clusters==0){
      throw std::logic_error("ERROR: number of clusters should not be zero\n");
  }
  if(m_num_dimensions==0){
      throw std::logic_error("ERROR: number of dimensions should not be zero\n");
  }
  for (j = 0; j < m_num_dimensions; j++) {
      sqsum += pow( ((*m_data)(k,j) - (*m_cluster_center)(i,j)) ,2) ;
  }
  return sqrt(sqsum);
}

double FCM::compute_membership_point(long i, long k){
    /*
     * i the cluster
     * k is the data point
    */
    //cout << __func__ <<"  num of cluster: "<<m_num_clusters<<endl;
    long j;
    double t, seg=0.0;
    double exp = 2 / (m_m - 1);
    double dik, djk;
    if(m_num_clusters==0){
        throw std::logic_error("ERROR: number of clusters should not be zero\n");
    }
    for (j = 0; j < m_num_clusters; j++) {
      dik = this->get_dist(i, k);
      djk = this->get_dist(j,k);
      if(djk==0){
          djk = DBL_MIN;
      }
      t = dik / djk;
      t = pow(t, exp);
      //cout << "cluster: " << i << "data: " << k << " - " << "t: "<<t<<endl;
      seg += t;
    }
    //cout << "seg: "<<seg << " u: "<<(1.0/seg)<<endl;
    return 1.0 / seg;
}


void FCM::set_data(MatrixXf *data){
    if(m_data!=nullptr){
        delete m_data;
    }
    if(data->rows()==0){
        throw std::logic_error("ERROR: seting empty data");
    }
    m_data = data;
    m_num_dimensions = m_data->cols();
}

void FCM::set_membership(MatrixXf *membership){
    if(m_data==0){
        throw std::logic_error("ERROR: the data should present before setting up the membership");
    }
    if(m_num_clusters==0){
        if(membership->cols() == 0){
            throw std::logic_error("ERROR: the number of clusters is 0 and the membership matrix is empty");
        }
        else{
            this->set_num_clusters(membership->cols());
        }
    }
    if(m_membership!=nullptr){
        delete m_membership;
    }
    m_membership = membership;
    if(m_membership->rows()==0){
        m_membership->resize(m_data->rows(), m_num_clusters);
    }
}

void FCM::init_membership(){
    long i, j;
    double mem;
    if(m_num_clusters == 0){
        throw std::logic_error("ERROR: the number of clusters is 0");
    }
    if(m_data==nullptr){
        throw std::logic_error("ERROR: the data should present before setting up the membership");
    }
    if(m_membership!=nullptr){
        delete m_membership;
    }
    m_membership = new MatrixXf;
    m_membership->resize(m_data->rows(), m_num_clusters);
    mem = 1.0 / m_num_clusters;
    for(j=0;j<m_num_clusters;j++){
        for(i=0;i<m_data->rows();i++){
            (*m_membership)(i,j) = mem;
        }
    }
}

void FCM::set_num_clusters(long num_clusters){
    m_num_clusters = num_clusters;
    if(m_cluster_center){
        delete m_cluster_center;
    }
    m_cluster_center = new MatrixXf;
    m_cluster_center->resize(m_num_clusters, m_num_dimensions);
}

MatrixXf * FCM::get_data(){
    return m_data;
}

MatrixXf * FCM::get_membership(){
    return m_membership;
}

MatrixXf * FCM::get_cluster_center(){
    return m_cluster_center;
}








