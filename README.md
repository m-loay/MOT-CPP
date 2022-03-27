# Multi Object Tracking

---

## Other Important Dependencies
* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 9.0
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* googletest [click here for installation instructions](https://github.com/google/googletest)
* gocavr 5.0 [click here for installation instructions](https://gcovr.com/en/stable/installation.html)

### Accuracy
The output RMSE px, py, vx, vy output coordinates [0.04852,0.062393,0.366935,0.400275].
Also the NIS for both lidar and radar is consistent as shown in the following figures:

For Lidar:

![lidar](media/highway--NIS_laser.png)


For Radar:

![radar](media/highway--NIS_Radar.png)

## Project results:

### UKF results:
For Lidar:

![lidarUKF](media/UKF_sample-laser-radar-measurement-data-1--NIS_laser.png)

For Radar:

![radarUKF](media/UKF_sample-laser-radar-measurement-data-1--NIS-Radar.png)

Output position Estimation VS Measurements VS GT:

![outputUKF](media/UKF_sample-laser-radar-measurement-data-1--Output-estimation.png)

Output velocity Estimation VS GT:

![velUKF](media/UKF_sample-laser-radar-measurement-data-1--Output-V.png)

### EkF results:
For Lidar:

![lidarEKF](media/EKF_sample-laser-radar-measurement-data-1--NIS_laser.png)

For Radar:

![radarEKF](media/EKF_sample-laser-radar-measurement-data-1--NIS-Radar.png)

Output position Estimation VS Measurements VS GT:

![outputEKF](media/EKF_estimation.png)

Output velocity Estimation VS GT:

![velEKF](media/EKF_sample-laser-radar-measurement-data-1--Output-V.png)

## Project Quality Work:

### Test Report Generated from google-test frame work:
![tc](media/test_cases.png)

### Coverage Reprot Report Generated from google-test frame work:
![cr](media/coverage.png)


gcovr -r . --exclude-throw-branches --html --html --html-details -o example-html-details.html

