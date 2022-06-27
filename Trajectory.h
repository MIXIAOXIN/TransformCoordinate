#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdint.h>
#include <vector>
#include <Eigen/Eigen>
#include "CoordinateSystemTrans.h"

namespace utility
{
	struct TrajPt
	{
		double _time_;
		double _x_;
		double _y_;
		double _z_;
		double _velocity_;
		TrajPt() {}
		TrajPt(double x, double y, double z)
		{
			_x_ = x;
			_y_ = y;
			_z_ = z;
		}
	};

	struct PosT
	{
		int seqNum;
		double gPSTime;
		double Northing;
		double Easting;
		double height;
		double latitude;
		double longitude;
		double speed;
		double roll;
		double pitch;
		double heading;
		std::string proName;
		int posQualtiy;
	};

	class Trajectory 
	{
	public:
		Trajectory();
		~Trajectory();

		void readTrajFile(const char* file_name);
		bool readPosTFile(const std::string filename, const float d_interval);
		bool readPosTFile_wuhan(const std::string filename, const float d_interval);
		bool readPosTFile_guangxi(const std::string filename);

		bool addTrajPts(TrajPt pt);
		//bool loadFromP2MTrajectory(P2M::BaseType::Traj3D::Ptr trajectory, const float traj_interval);

		void writeTraj_XYZT(const std::string file_name);
		void writeTrajToLas(const char* file_name);
		inline std::vector<TrajPt> getTrajPts()
		{
			return _traj_pts_;
		}

		bool getTrajPtsInTimeInterval(const std::vector<double>& time_interval,
			std::vector<TrajPt>& traj_pts);

		static void writeTraj(const std::string& file_name,
			const std::vector<TrajPt>& traj_pt,
			const int start_id,
			const int end_id);



		inline
		float calculateDistance2D(const Eigen::Vector4d &pt1, const Eigen::Vector4d &pt2) {
			float xx, yy, dist;
			xx = pt1[0] - pt2[0];
			yy = pt1[1] - pt2[1];
			dist = sqrtf(xx*xx + yy * yy);
			return dist;
		}

		inline
			float calculateDistance2D(const TrajPt &pt1, const TrajPt &pt2) {
			float xx, yy, dist;
			xx = pt1._x_- pt2._x_;
			yy = pt1._y_ - pt2._y_;
			dist = sqrtf(xx*xx + yy * yy);
			return dist;
		}
	private:
		std::vector<TrajPt> _traj_pts_;

	private:
		/**
		* \brief: binary search for the given time stamp,
		* \return: trajectory id in the trajectory vector.
		*/
		int binarySearchForGivenTime(double time);
	};

} // namespace utility

#endif



