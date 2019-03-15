#ifndef LANCZOS_HPP
#define LANCZOS_HPP

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <lapacke.h>
#include <cblas.h>
#include <random>
#include <chrono>

template <typename MatrixType>
class LanczosSolver{
public:
	enum class EigenvalueType{};

	typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::Index Index;
    typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> DenseMatrix;
    typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, 1> ScalarVector;

    LanczosSolver() = delete;
    LanczosSolver(const DenseMatrix& a, Index num_singular_values);
    LanczosSolver(const LanczosSolver& other) = delete;
    LanczosSolver(LancosSolver&& other);
    ~LanczosSolver();

    void compute(Index num_eigen_values);
    void compute(const DenseMatrix& a, Index num_singular_values);

    const DenseMatrix& matrixU();
    const DenseMatrix& matrixV();
    const ScalarVector& singularValues();

private:

	ScalarVector modifedGramSchmidt(const DenseMatrix& mat, Index num_cols, const ScalarVector& v){
		LanczosSolver<MatrixType>::ScalarVector w = v;

		if(num_cols == 0){
			return w / w.norm();
		}

		for(Index i = 0; i < num_cols; ++i){
			float h = w.dot(mat.col(i));
			w = w - h * mat.col(i);
			if(w.norm() < std::numeric_limits<MatrixType::Scalar>::epsilon()){
				break;
			}
		}

		if(w.norm() < std::numeric_limits<MatrixType::Scalar>::epsilon()){
			return LanczosSolver<MatrixType>::ScalarVector::Zero(w.size());
		}
		else{
			return w/w.norm();
		}
	}

	enum class Mode{};

	void implicitRestartArnoldi(Mode mode, const DenseMatrix& a, EigenvalueType type, Index num_ev, Scalar tolerance, DenseMatrix& krylov_, Index max_iterations,
		){

	}

	void computeEigenVectors(){

	}


private:
	DenseMatrix krylov_;
	std::mt19937 rng_;
}



#endif