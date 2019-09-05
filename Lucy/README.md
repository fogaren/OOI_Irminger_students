
## This document explains all the code used in Lucy Wanzer's thesis. See Master_thesis_Wanzer for an overview of all the scripts and high level comments. See each individual script for more high level comments. This readme file will give detailed comments for each script. 

### Script 1: Irminger005_CTD
#### Irminger005_CTD loads in all the CTD data from the Irminger 005 cruise. Lines 3 and 4 identify where the excel sheet is stored on one's computer and then identifies the file name. Line 6 pulls out which casts we worked with in this script. We were looking at casts that were deeper than 1500 m. The for loop (lines 7-30) assigns the correct parameter to each column for every cast that we picked out. This is done by  telling MATLAB which parameter is in each column and including all rows. 

#### The majority of the rest of the script (until line 164) is just playing with the data and plotting various parameters of various casts together so that I could identify the differences between the different casts. The section that is line 164 to 246 identifies the casts with the shallowest winter ventilation depths and the casts with the deepest winter ventilation depths, using both AOU and density. Lines 212 to 246 plots the shallower and deeper CTD casts next to each other on AOU and density profiles. Lines 268 to 300 plot the locations of the OOI moorings and the CTD casts. The color of the CTD casts indicates the depth of winter ventilation. 

### Script 2: load_HYPM
#### load_HYPM loads in 4 years of wire-following profile data from the profiler mooring. Lines 3-156 load in the data from each year. The same steps are repeated for each year, so I will use Year 1 as an example. Line 5 identifies the filename and also uses ncdisp function to display all the nitty gritty information about the data file. When you run this section, all the information about the particular dataset including names of individual parameters, how large the file is etc. will show up in your command window. This is important because the names of parameters often change between years. For example, in Year 1, the name for salinity is much longer than in Year 4. Using ncdisp will insure that you're using the correct name for each year. This section (lines 4-156) extracts all the data and assigns variable names to each type of data. 

#### Lines 158-186 use functions from the Gibbs seawater toolbox to make conversions between parameters to get the correct variables which will later be used to grid in various ways. Lines 189-260 grids the data onto depth and temperature grids. This section starts by defining the interval across which the data will be gridded. For example in the depth grid, the profiles will be compiled from 150 to 2600 m, and will be gridded every 5 meters. The next section takes the mean of the up and down profiles for all the variables (scivars) for both depth and temperature. 

#### Lines 279 to 361 assigns the scivars to new variable names that identify them by either depth or temperature grid and what they are (e.g. salinity, temp etc). Using the squeeze function, the profiles are put into the correct dimensions to work with. Lines 364-382 use Gibbs seawater toolbox functions to calculate oxygen saturation.

#### Lines 383-420 creates 4 figures (one figure for each year of data) and each figure contains initial timeseries plots for all the parameters listed. 

#### Lines 422-431 merge all four years of data into a single dataset. In the following section, this full dataset is plotted in order to visualize the full dataset all together. In the final section (433-497), X and Y are assigned as time and depth which is then used to plot all parameters for the full four years of data. At the top of this section are adjustable parameters like min and max depth and the color map -- by changing variables here, they will change in the plots below. 

### Script 3: wfp_Irminger_winklercalibration_allyrs
#### This script was mainly written by Hilary so I will mostly give a general overview. The first thing that this script does is run another script called "loadWinklerIrmingerAllYears" which loads in all the discrete oxygen data needed to make the initial gain corrections to the dissolved oxygen data from the profiler mooring. To run this script you will need to have 4 excel sheets in your path which contain all the discrete oxygen data. All four excel sheets are contained in the folder "CruiseData_Yrs1to4" which is in the OOI_Irminger_students folder. Each excel sheet is named by the year the data is coming from and the cruise the data was collected from, for example, AR30 and KN221 denote the research vessels used and the specific cruise number. 

#### For each year, after reading in the excel data, all the columns are assigned to Matlab variables and then the standard deviation is calculated for each year and the mean standard deviation is used. All the samples that are duplicates are also identified which then allows us to calculate the error. 

#### Back to Script 3... lines 6-9 identify the parameters that are decided upon for comparing the discrete oxygen samples to the oxygen data taken from the profiler mooring. For example, line 7 describes the depth_window which is the limit of the distance that we compare profiles and discrete samples. The limitations in these four lines can be adjusted depending on preference. 

#### In lines 11-27, casts are identified from each year that are in the OOI site region. I adjusted this list by going back through all the spreadsheets and finding which casts were in the OOI region and then adding them to this list. Lines 61-63 assign variables to the gain correction number and the gain standard deviation for each year. Make sure to go to your workspace and look at these numbers when running the script to make sure they make sense. The rest of the script makes plots of what the gain correction looks like for each cast in each year period. 

### Script 4: HYPM_O2_gaincorr
#### Lines 3-18 creates new variables for each year (denoted with "_corr") of gain-corrected oxygen values by multiplying the oxygen values in the depth grid, temperature grid (isotherms), and unguided data. 

#### Lines 22-25 calculate O2 saturation with the new gain corrected data



### Script 5: wfp_deepIsotherm_driftCorrection
#### This script does all of the calculations for how the oxygen data drifts overtime. The instruments on the profiler mooring become less accurate overtime and so this script picks out 5 stable deep isotherms (these are shown on line 6) and looks at how the oxygen data for these isotherms changes throughout the year. Theoretically 
