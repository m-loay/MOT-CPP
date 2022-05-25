/*
 * kf_lib.h
 *
 *  Created on: Apr 18, 2019
 *      Author: mody
 */
/** @file kf_lib.h
 *  @ingroup kf
 *  @brief kf class
 */

/**
 *  @addtogroup kf
 *  @{
 */

#ifndef KALMAN_FILER_H_
#define KALMAN_FILER_H_
#include <Eigen/Dense>
#include <vector>
#include "KalmanConfig.h"

class kalmanFilter
{
public:
	kalmanFilter(int num_states)
	{
        //save number of states
        m_numState = num_states;

        //create a 4D state vector, we don't know yet the values of the x state
        x = Eigen::VectorXd(m_numState);
        x.fill(0.0);

        //state covariance matrix P
        P = Eigen::MatrixXd::Identity(m_numState, m_numState);
        P.diagonal().fill(0.0);
	}

	kalmanFilter(int num_states,Eigen::VectorXd state):kalmanFilter(num_states)
	{
		x = state;
	}

	kalmanFilter(int num_states,
			     Eigen::VectorXd state,
				 Eigen::MatrixXd Covar):kalmanFilter(num_states, state)
	{
        P = Covar;
	}

	Eigen::VectorXd getMean(void)const
	{
		return x;
	}

	Eigen::MatrixXd getCovariance(void)const
	{
		return P;
	}

	/**
	 * @brief predict
	 * Perform the Prediction step in kalman filter.
	 *
	 * @param[in,out] x
	 * The mean of stat vector  {VectorXd &}.
	 *
	 * @param[in,out] P
	 * The mean of stat vector  {VectorXd &}.
	 *
	 * @param[in] Q
	 * the process noise matrix  {MatrixXd}.
	 *
	 * @param[in] g
	 * calculates the mean state vector based dynamic model{function}.
	 *
	 * @param[in] g_prime
	 * calculates the mean state vector based dynamic model{function}.
	 *
	 * @param[in] p_args
	 *  Extra arguments {const void *}.
	 *
	 */
	void predict(const Eigen::Ref<const Eigen::MatrixXd> &Q,
			     std::function<Eigen::VectorXd(const Eigen::Ref<const Eigen::VectorXd> &, const void *)>g,
			     std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::VectorXd> &, const void *)>g_prime,
			     const void *p_args = NULL)
	{

		x = g(x, p_args);
		Eigen::MatrixXd G = g_prime(x, p_args);
		P = (G* P * G.transpose()) + Q;
	}

	/**
	 * @brief update
	 * Perform the update step.
	 *
	 * @param[in,out] x
	 * The mean of stat vector  {VectorXd &}.
	 *
	 * @param[in,out] P
	 * The mean of stat vector  {VectorXd &}.
	 *
	 * @param[in] Y
	 * the innovation vector  {VectorXd}.
	 *
	 * @param[in] H
	 * the measurement matrix  {MatrixXd}.
	 *
	 * @param[in] k
	 * the kalman gain matrix  {MatrixXd}.
	 *
	 */
	void update(const Eigen::Ref<const Eigen::VectorXd> &Y,
			    const Eigen::Ref<const Eigen::MatrixXd> &H,
			    const Eigen::Ref<const Eigen::MatrixXd> &R)
	{

		CalculateKalmanGain(H,R);
		x = x + (K * Y);
		Eigen::MatrixXd I(Eigen::MatrixXd::Identity(x.size(), x.size()));
		P = (I - K * H) * P;
	}

	/**
	 * @brief update
	 * Perform the update step.
	 *
	 * @param[in,out] x
	 * The mean of stat vector  {VectorXd &}.
	 *
	 * @param[in,out] P
	 * The mean of stat vector  {VectorXd &}.
	 *
	 * @param[in] Y
	 * the innovation vector  {VectorXd}.
	 *
	 * @param[in] H
	 * the measurement matrix  {MatrixXd}.
	 *
	 * @param[in] k
	 * the kalman gain matrix  {MatrixXd}.
	 *
	 */
	void update(const Eigen::Ref<const Eigen::VectorXd> &Y,
			    const Eigen::Ref<const Eigen::MatrixXd> &H,
			    const Eigen::Ref<const Eigen::MatrixXd> &R,
				std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd> &,
						                      const Eigen::Ref<const Eigen::MatrixXd> &)>kFun)
	{

		kFun(H,R);
		x = x + (K * Y);
		Eigen::MatrixXd I(Eigen::MatrixXd::Identity(x.size(), x.size()));
		P = (I - K * H) * P;
	}
private:
	/**
	 * @brief CalculateKalmanGain
	 * It calculates kalman gain.
	 *
	 * @param[in] P
	 * The mean of stat vector  {VectorXd &}.
	 *
	 * @param[in] H
	 * the measurement matrix  {MatrixXd}.
	 *
	 * @param[in] R
	 * the noise covariance measurement matrix  {MatrixXd}.
	 *
	 * @return K
	 * the kalman gain matrix  {MatrixXd}.
	 *
	 */
	void CalculateKalmanGain(const Eigen::Ref<const Eigen::MatrixXd> &H,
			                 const Eigen::Ref<const Eigen::MatrixXd> &R)
	{
		Eigen::MatrixXd Ht(H.transpose());
		Eigen::MatrixXd S((H * P * Ht) + R);
		K = (P *Ht * S.inverse());
	}

    Eigen::VectorXd x;	// object state
    Eigen::MatrixXd P;	// object covariance matrix
    Eigen::MatrixXd K;
    int m_numState;
};

#endif /* KALMAN_FILER_H_ */
/**
 *  @}
 */
