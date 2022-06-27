#include "Trajectory.h"

namespace utility
{

	Trajectory::Trajectory()
	{
		_traj_pts_.resize(0);
	}

	Trajectory::~Trajectory()
	{
	}

	void Trajectory::readTrajFile(const char *file_name)
	{
		FILE *fp = fopen(file_name, "r");
		if (fp == NULL)
		{
			printf("open file failed��%s", file_name);
			return;
		}
		while (!feof(fp))
		{
			TrajPt tp;
			fscanf(fp, "%lf,%lf,%lf\n", &tp._x_, &tp._y_, &tp._z_
				   //, &tp._time_
			);
			// std::cout << tp._x_ << ',' << tp._y_ << ',' << tp._z_ << ',' << tp._time_ << std::endl;
			_traj_pts_.push_back(tp);
		}
		fclose(fp);
	}

	bool Trajectory::readPosTFile(const std::string filename, const float d_interval)
	{
		/*****************interval for adjacent trajectory points*****************/
		float interval = d_interval;
		/**********************************/

		std::ifstream inf;
		inf.open(filename);
		if (!inf.is_open())
		{
			std::cout << "failed to open the .PosT file!" << filename << std::endl;
			return false;
		}

		std::string line;
		for (size_t i = 0; i < 25; i++)
		{
			std::getline(inf, line); // PosT
		}

		int traj_pts_count = 0;
		Eigen::Vector4d currentPt, lastPt;
		while (!inf.eof())
		{
			if (inf.fail())
				break;

			float dist;
			PosT post;
			inf >> post.seqNum >> post.gPSTime >> post.Northing >> post.Easting >> post.height >> post.latitude >> post.longitude >> post.speed >> post.roll >> post.pitch >> post.heading >> post.proName >> post.posQualtiy;
			currentPt[0] = post.Easting;
			currentPt[1] = post.Northing;
			currentPt[2] = post.height;
			currentPt[3] = post.gPSTime;
			TrajPt tp;
			tp._x_ = currentPt[0];
			tp._y_ = currentPt[1];
			tp._z_ = currentPt[2];
			tp._time_ = currentPt[3];
			if (traj_pts_count == 0 || post.seqNum == 0)
			{
				_traj_pts_.push_back(tp);
				lastPt = currentPt;
				++traj_pts_count;
			}
			else
			{
				dist = calculateDistance2D(lastPt, currentPt);
				if (dist > interval)
				{
					_traj_pts_.push_back(tp);
					lastPt = currentPt;
					++traj_pts_count;
				}
			}
		}
		inf.close();
	}

	bool Trajectory::readPosTFile_wuhan(const std::string filename, const float d_interval)
	{
		/*****************interval for adjacent trajectory points*****************/
		float interval = d_interval;
		/**********************************/

		std::ifstream inf;
		inf.open(filename);
		if (!inf.is_open())
		{
			std::cout << "failed to open the .PosT_wuhan file!" << filename << std::endl;
			return false;
		}

		std::string line;
		for (size_t i = 0; i < 1; i++)
		{
			std::getline(inf, line); // PosT file instructons: the head lines for file introduction
		}

		Eigen::Vector4d currentPt, lastPt;
		int count = 0;
		while (!inf.eof())
		{
			if (inf.fail())
				break;

			float dist;
			short it_year, it_month, it_day, it_hour, it_minute, it_second, it_mmsecond;
			PosT post;
			inf >> post.seqNum >> post.gPSTime >> it_year >> it_month >> it_day >> it_hour >> it_minute >> it_second >> it_mmsecond >> post.Easting >> post.Northing >> post.height >> post.heading >> post.pitch >> post.roll;
			currentPt[0] = post.Easting;
			currentPt[1] = post.Northing;
			currentPt[2] = post.height;
			currentPt[3] = post.gPSTime;
			TrajPt tp;
			tp._x_ = currentPt[0];
			tp._y_ = currentPt[1];
			tp._z_ = currentPt[2];
			tp._time_ = currentPt[3];
			if (post.seqNum == 0 || _traj_pts_.size() < 1)
			{
				_traj_pts_.push_back(tp);
				lastPt = currentPt;
				++count;
			}
			else
			{
				dist = calculateDistance2D(lastPt, currentPt);
				if (dist > interval)
				{
					_traj_pts_.push_back(tp);
					lastPt = currentPt;
					++count;
				}
			}
		}
		inf.close();
		std::cout << "Selected trajectory point size: " << count << std::endl;
	}

	bool Trajectory::readPosTFile_guangxi(std::string filename)
	{
		std::ifstream inf;
		inf.open(filename);
		if (!inf.is_open())
		{
			std::cout << "failed to open the .PosT_guangxi file!" << filename << std::endl;
			return false;
		}

		std::string line;
		for (size_t i = 0; i < 25; i++)
		{
			std::getline(inf, line); // PosT file instructons: the head lines for file introduction
		}

		Eigen::Vector4d currentPt, lastPt;
		int count = 0;
		while (!inf.eof())
		{
			if (inf.fail())
				break;

			short it_year, it_month, it_day, it_hour, it_minute, it_second, it_mmsecond;
			float veast, vnorth, vup, sdnorth, sdeast, sdheight, henghd, pithsd, rollsd;
			PosT post;
			inf >> post.gPSTime >> post.latitude >> post.longitude >> post.height >>
				post.heading >> post.pitch >> post.roll >>
				vnorth >> veast >> vup >> sdnorth >> sdeast >> sdheight >> henghd >> pithsd >> rollsd >> post.seqNum;

			// coordinate transformation:
			TrajPt tp;
			WGS84ToGauss3(post.latitude, post.longitude, post.height, 111.0, post.Easting, post.Northing, post.height);
			//LongtitudeLatitude2UTM(post.latitude, post.longitude, post.Easting, post.Northing);

			tp._x_ = post.Easting;
			tp._y_ = post.Northing;
			tp._z_ = post.height;
			tp._time_ = post.gPSTime;
			_traj_pts_.push_back(tp);
			++count;
		}
		inf.close();
		std::cout << "Selected trajectory point size: " << count << std::endl;
	}

	void Trajectory::writeTraj_XYZT(const std::string file_name)
	{
		std::ofstream ofs;
		ofs.open(file_name, std::ios::out | std::ios::binary);
		if (!ofs.is_open())
		{
			printf("write file failed:%s", file_name);
			return;
		}

		try
		{
			for (int i = 0; i < _traj_pts_.size(); i++)
			{
				ofs << std::fixed << std::setprecision(3) << _traj_pts_[i]._x_ << ','
					<< _traj_pts_[i]._y_ << ','
					<< _traj_pts_[i]._z_ << ','
					<< _traj_pts_[i]._time_ << std::endl;
			}
			ofs.flush();
			ofs.close();
		}
		catch (std::exception *e)
		{
			ofs.close();
		}
	}

	void Trajectory::writeTraj(const std::string &file_name,
							   const std::vector<TrajPt> &traj_pt,
							   const int start_id,
							   const int end_id)
	{
		std::ofstream ofs;
		ofs.open(file_name, std::ios::out | std::ios::binary);
		if (!ofs.is_open())
		{
			printf("write file failed:%s", file_name);
			return;
		}

		try
		{
			for (int i = start_id; i <= end_id; i++)
			{
				ofs << std::fixed << std::setprecision(3) << traj_pt[i]._x_ << ','
					<< traj_pt[i]._y_ << ','
					<< traj_pt[i]._z_ << std::endl;
			}
			ofs.flush();
			ofs.close();
		}
		catch (std::exception *e)
		{
			ofs.close();
		}
	}

	bool Trajectory::getTrajPtsInTimeInterval(const std::vector<double> &time_interval,
											  std::vector<TrajPt> &traj_pts)
	{
		bool re_flag = true;
		std::vector<TrajPt>().swap(traj_pts);
		int nums_traj_points = _traj_pts_.size();
		int start_index = binarySearchForGivenTime(time_interval[0]);
		int end_index = binarySearchForGivenTime(time_interval[1]);
		if (start_index >= end_index)
		{
			std::cout << "time interval begin end wrong !!!!" << std::endl;
			re_flag = false;
			return re_flag;
		}

		// Eigen::Vector4d first_pt;
		// Eigen::Vector4d current_pt;
		// first_pt[0] = (_traj_pts_.begin() + start_index)->_x_;
		// first_pt[1] = (_traj_pts_.begin() + start_index)->_y_;
		// first_pt[2] = (_traj_pts_.begin() + start_index)->_z_;
		// first_pt[3] = (_traj_pts_.begin() + start_index)->_time_;
		// traj_pts.push_back(_traj_pts_[start_index]); // insert the first trajectory point

		// float max_distance = 0.0f;
		// for (int i = start_index + 1; i < end_index; ++i)
		//{
		//	current_pt[0] = (_traj_pts_.begin() + i)->_x_;
		//	current_pt[1] = (_traj_pts_.begin() + i)->_y_;
		//	current_pt[2] = (_traj_pts_.begin() + i)->_z_;
		//	current_pt[3] = (_traj_pts_.begin() + i)->_time_;
		//	float current_distance = calculateDistance2D(first_pt, current_pt);

		//	// insert the trajectory point in increasing distance
		//	if (current_distance > max_distance)
		//	{
		//		max_distance = current_distance;
		//		traj_pts.push_back(_traj_pts_[i]);
		//	}
		//}

		traj_pts.assign(_traj_pts_.begin() + start_index, _traj_pts_.begin() + end_index);
#ifdef debug_trajectory_selection
		std::string traj_selected_filename = "trajectory_pts_selected.txt";
		this->writeTraj(traj_selected_filename, _traj_pts_, start_index, end_index);
#endif
		return re_flag;
	}

	int Trajectory::binarySearchForGivenTime(double time)
	{
		int min = 0;
		int max = _traj_pts_.size() - 1;
		int mid = min;

		while (min < max - 1)
		{
			mid = (min + max) / 2;
			if (_traj_pts_[mid]._time_ > time)
			{
				max = mid;
			}
			else
			{
				if (_traj_pts_[mid]._time_ == time)
				{
					break;
				}
				min = mid;
			}
		}
		return mid;
	}

	bool Trajectory::addTrajPts(TrajPt pt)
	{
		_traj_pts_.push_back(pt);

		return true;
	}

} // namespace utility