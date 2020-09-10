// Supernova Dataset Manager

#include <iostream>
#include <string>
#include <vector>
#include <iterator>
#include <math.h>
#include <fstream>
#include <algorithm>
#include <experimental/filesystem>

using namespace std;

class sn {
    static bool lums;      // false for absolute magnitudes, true for luminosities (user-specified)
    float lumdist;        // The luminosity distance in parsecs
    vector<float> t_peak;  // the time of peak luminosity in each selected band
  public:
    static int sn_count;
    sn(bool choice) {sn_count++; lums = choice;};  // constructor, only called first time
    string sn_name;
    string sn_type;
    sn(string & name,string & type) {sn_count++; sn_name = name; sn_type = type;};   // other constructor, called the rest of times
    // the following 2D vectors will store our data to be output, the first index will be the # of photometric bands
    // (see below) while the 2nd index will be for the different observations in that band
    vector<vector<float>> time;  // the time of observations (read in as MJD and later converted to days relative
                                 // to peak luminosity
    vector<vector<double>> lum;   // originally read-in as apparent magnitude and later converted to absolute
                                   // magnitude or luminosities in (erg/s) (see the private variable lums)
    void find_lumdist(ifstream&);  // sets the private lumdist variable from reading the file
    void build_time_and_lum(ifstream&, vector<string>&); // build the time and luminosity 2D vectors.
    void convert_mags();   // converts lum 2D vector from apparent magnitude to absolute magnitude or luminosity in erg/s
    void find_peak_times();  // finds the times of peak luminosity
   						// for each photometric band, building private vector t_peak 
    void offset_time_by_peak();    // converts from the MJD input to days relative to t_peak
};

bool sn::lums; // Just to initialize, will be set based on user input

int sn::sn_count = -1; // increments everytime an object of class sn is created
		       // begins at -1 because there will be at least one call
		       // to the constructor that sets lums.

void sn::find_lumdist(ifstream & stream){

    string line;
    size_t pos;
  
    bool found = false;
    string result = "";
    while(getline(stream, line)){
        if (line.find("\"lumdist\"") != string::npos){
            while(getline(stream, line)){
                if (line.find("\"value\"") != string::npos){
                    pos = line.find("value") + 5; 
          	       
                    while (!isdigit(line[pos])) {  // move pos until we find a digit
   	                ++pos;
                    }
                    while (line[pos] != '"'){  // the lumdist's value end is signaled by closing quotes
                        result += line[pos];
                        ++pos;
  	            }
  	            found = true;
  	            break; // we have our answer, just need to convert to a float (and Mpc to pc)  
                }
            }	     
        }
        if (found) break;
    } 
    lumdist = stof(result) * pow(10,6); 
    return;
}

void sn::build_time_and_lum(ifstream & stream, vector<string> & req_sn_bands){
  
    string line;
    size_t pos;
  
    string t, b, m;   // temporary strings to store single values of time, band, and magnitude
    t = "";
    m = "";
  
    vector<string>::iterator it;
    int i;
    int j [req_sn_bands.size()];
  
    vector<float> tempf;
    vector<double> tempd;
   
    // this loop necessary to avoid seg faults;
    // In the event that there is no matching data
    // for one of the requested bands i, the vectors
    // time[i] and lum[i] are single length.
    
    for (int i = 0; i < req_sn_bands.size(); ++i) { //create initial space in the 2D vectors
        j[i] = 0;
        time.push_back(tempf);
        time[i].push_back(0.f);
        lum.push_back(tempd);
        lum[i].push_back(0.);
    }
  
    while(getline(stream, line)){
        if (line.find("\"photometry\"") != string::npos){
            while(getline(stream, line)){
                if (line.find("]") != string::npos) break;  // if we reach the end of the 'block'
  	                                                    // titled "photometry" in the file, break
                // The following conditional blocks will cycle in order until the end of the block. So
  	        // 2D vector extensions and the filling of their elements happens after magnitude is read
  	        // for each observation.
  	        if (line.find("\"time\"") != string::npos){
                    pos = line.find("time") + 7;  // length of the word 'time' + # of chars before value starts
                    while (line[pos] != '"'){
  	                t += line[pos];	 
  	                ++pos; 
                    }
                }
  	        else if (line.find("\"band\"") != string::npos){
  	            pos = line.find("band") + 7;  // length of the word 'band' + # of chars before value starts
                    while (line[pos] != '"'){
  	                b += line[pos];
  	                ++pos;
                    }	  
  	            // if band does not match continue to next observation
  	            it = find(req_sn_bands.begin(), req_sn_bands.end(),b);
                    if (it == req_sn_bands.end()) {
          	        t = "";
          	        b = "";    
          	        continue; 
          	    }
          	    else i = it - req_sn_bands.begin();   // the row index to edit for time and luminosity
          	}
                else if (line.find("\"magnitude\"") != string::npos){
  	            if (t == "") continue;  // skip magnitude if the last read-in band does not match
  	            pos = line.find("magnitude") + 12;  // length of the word 'magnitude' + # of chars before value starts
                    while (line[pos] != '"'){
                        m += line[pos];
  	                ++pos;
  	            }
  	           // fill appropriate 2D vector elements (allocating space when needed) and re-initialize temporary strings
          	   while (j[i] > time[i].size()-1) {time[i].push_back(0.f); lum[i].push_back(0.);}
         	   time[i][j[i]] = stof(t);
         	   lum[i][j[i]] = stod(m);
                   ++j[i];  // increment column index for band (row) i 	  
         	   t = "";
         	   b = "";
         	   m = "";
                }	
            }	      
        }
    }  
}

void sn::convert_mags(){

    for (int i = 0; i < lum.size(); ++i){
        for (int j = 0; j < lum[i].size(); ++j){
            lum[i][j] -= 5 * (log10(lumdist) - 1);  // absolute magnitude
            if (lums){
                lum[i][j] /= -2.5;
                lum[i][j] = pow(10, lum[i][j]);
                lum[i][j] *= 3.0128 * pow(10,35);   // luminosity in erg/s
            }
        }
    }
}

void sn::find_peak_times(){
	
    for (int i = 0; i < time.size(); ++i){	// loop over photometric bands  
        double max_lum = 0.;
        int index = 0;
        t_peak.push_back(0.f);
        for (int j = 0; j < time[i].size(); ++j){      // loop over observations
            if (abs(lum[i][j]) > abs(max_lum)){
                max_lum = lum[i][j];
                index = j;
            }
        }
        t_peak[i] = time[i][index];
    }  
}

void sn::offset_time_by_peak(){
  
    for (int i = 0; i < time.size(); ++i){
      for (int j = 0; j < time[i].size(); ++j){
        time[i][j] -= t_peak[i]; 
      }
    }
}

string get_name(ifstream & stream){

    string line;
    size_t pos;
  
    string result = "";
    while(getline(stream, line)) {
        if (line.find("\"name\"") != string::npos){
          pos = line.find("name") + 4;  // move pos to end of the word "name"
          while (!isalpha(line[pos])){  // move pos until we find an alpha char (starts the name)
              ++pos;       
          }
          while (line[pos]!='"'){  // the name's end is signaled by the closing quotes
              result += line[pos];
              ++pos;
          }
          return result;
        }
    }
    return result;  // return empty string if the file does not have the name field
}

string get_type(ifstream & stream, vector<string> & req_sn_types){

    // Function to retrieve and return the SN type from the file 
    string line;
    size_t pos;

    string result = "";

    while(getline(stream, line)) {
        if (line.find("\"claimedtype\"") != string::npos){
            while (getline(stream, line)) {
                if (line.find("]") != string::npos) break;  // if we reach the end of the 'block'
                                                           // titled "claimedtype" in the file, break
                if (line.find("\"value\"") != string::npos){
                    pos = line.find("value") + 8;          // the length of the word "value" + # of chars before
    	  	            				 // value starts
                    while (line[pos] != '"'){           // the type's end is signaled by closing quotes
                        if (line[pos] != ' ') result += line[pos]; // want no spaces to go into result
                            ++pos;
                        }
                        if (find(req_sn_types.begin(), req_sn_types.end(), result) == req_sn_types.end()) {
                            result = "";
                            continue;  // no match found for string s, check other values in file
    	                }  
                        else return result;
                    }
                }
                return result; // return empty string if no match found
            }
        }
    return result; // return empty string if the file does not have the claimedtype field
}

void fill_req(vector<string> & req, vector<string> & choices, const string & parameter){

    // Function to fill the arrays of requested types and bands (req_sn_types and req_sn_bands)
    // from user input.	
    vector<string>::iterator it;	
    while(req.size() < choices.size()){
        string s;
	cout << "Enter a " << parameter << " from the list, otherwise enter 'done':" << endl;
	cin >> s;
        if (s == "done") break;                                                                               
        it = find(choices.begin(), choices.end(), s);                                                         
        if (it == choices.end()){                                                                            
	    cout << "You entered an invalid " << parameter << ". Try again." << endl;                               
	    continue;                                                                                         
	}                                                                                                     
        it = find(req.begin(), req.end(), s);
	if (it != req.end()){
            cout << "You've already entered this " << parameter << ". Try again." << endl;
	    continue; 
	}
	req.push_back(s);
	cout << "spectral type " << s << " added to list of requested SN types" << endl;	
    }
    cout << "Your list of requested " << parameter << "s (";
    for (int i = 0; i < req.size() - 1; ++i){
        cout << req[i] << ", ";
    }
    cout << req[req.size() - 1] << ") has been stored." << "\n" << endl;
}  

int main()
{

    cout << "Welcome to the Supernova Data Manager!" << endl;
    cout << "The purpose of this program is to organize data" << endl;
    cout << "according to the user's wishes and output it for" << endl;
    cout << "subsequent analysis (plotting, etc.)" << "\n" << endl; 
    cout << "We will soon read data from input files" << endl;
    cout << "(with .jso extensions, available for download" << endl;
    cout << "from The Open Supernova Catalog (https://sne.space))" << endl;
    cout << "but first you will be asked questions to help us" << endl;
    cout << "pre-filter the data to be read" << "\n" << endl;
  
    vector<string> sn_types = {"Ia", "Ib", "Ic", "II", "IIn", "IIP", "IIPec"}; // allowed spectral types
    vector<string> req_sn_types;  // stores the requested SN spectral types
  
    vector<string> sn_bands = {"U","B","V","R","I"}; // allowed photometric bands
    // UV, Blue, Visible, Red, IR. Bands longer than one character do exist. 
    // Type is vector<string> for easier extension. Only 1 character bands included
    // because these happen to be the most commonly used.
    vector<string> req_sn_bands;  // stores the requested SN photometric bands
  
    vector<string>::iterator it;
  
    cout << "Select from the following list all the spectral types you wish to extract data for:" << endl;
    cout << "Ia,Ib,Ic,II,IIn,IIP,IIPec" << "\n" << endl;
    cout << "Supernova files which do not match any of your selected types" << endl;
    cout << "will not have their data stored." << endl;

    fill_req(req_sn_types, sn_types, "SN type");

    cout << "Select from the following list all the photometric bands you wish to extract data for:" << endl;
    cout << "U,B,V,R,I" << "\n" << endl;
    cout << "Supernova data which were not observed in any of your selected bands" << endl;
    cout << "will not be stored. The data that is stored will separately stored based" << endl;
    cout << "on the appropriate band. We will finally output, files which store the data" << endl;
    cout << "for each SN type and observed photometric band." << endl;
  
    fill_req(req_sn_bands, sn_bands, "spectral type");
/*    while(req_sn_bands.size() < sn_bands.size()){
        string s;
        cout << "Enter a photometric band from the list, otherwise enter 'done':" << endl;
        cin >> s;
        if (s == "done") break;
        it = find(sn_bands.begin(), sn_bands.end(), s);
        if (it == sn_bands.end()) {
            cout << "You entered an invalid photometric band. Try again." << endl;
            continue;
        }
        it = find(req_sn_bands.begin(), req_sn_bands.end(),s);
        if (it != req_sn_bands.end()) {
            cout << "You've already entered this photometric band. Try again." << endl;
            continue;
        }
        req_sn_bands.push_back(s);
        cout << "photometric band " << s << " added to list of requested bands" << endl;
    }
  
    cout << "Your list of requested SN bands (";
    for (int i = 0; i < req_sn_bands.size() - 1; ++i){
        cout << req_sn_bands[i] << ", ";
    }
    cout << req_sn_bands[req_sn_bands.size() - 1] << ") has been stored." << "\n" << endl;
*/

    cout << "Finally, before we begin reading data, would you like the apparent magnitudes" << endl;
    cout << "stored in the files to be converted to absolute magnitudes or luminosities" << endl;
    cout << "(in erg/s)? (Enter 'lums' for luminosities and anything else for abs. magnitudes):" << endl;
  
    string s;
    cin >> s;
  
    if (s == "lums") sn supernova(true);   // set the private bool var lums for the class sn
    else sn supernova(false);    
  
    vector<sn> supernovae;    // a vector to store objects of class sn
    string filename;

    for (auto& p: experimental::filesystem::directory_iterator(".")){
        filename = p.path();
        // check extension. continue if not .json
        if (filename.substr(filename.size() - 5, filename.size() - 1) != ".json") continue;
        ifstream snfile(filename);   // open file
        cout << "Reading file: " << filename << " ..." << endl;
        string name = get_name(snfile);  // get supernova name
        string type = get_type(snfile, req_sn_types); // get supernova type
        // close file and continue if SN type is not in list of requested types
        if (type == "") {
            cout << "None of the claimed types for this superova match any of your requested types. Closing file." << endl;
            snfile.close(); 
            continue;
        }
        cout << "Match found between reqested spectral types and this SN's type: " << type << endl;
        sn supernova(name, type); // constructor
        supernova.find_lumdist(snfile);
        supernova.build_time_and_lum(snfile, req_sn_bands);
        cout << "Closing file: " << filename << endl;
        snfile.close();
        supernova.convert_mags();
        supernova.find_peak_times();
        supernova.offset_time_by_peak();
        supernovae.push_back(supernova);
    }

    cout << "You have imported data for " << sn::sn_count << " supernova(e)" << endl;
    cout << "Now the program will output data to filenames based on your requested" << endl;
    cout << "spectral types and the names of their matching supernovae." << "\n" << endl;

    string dir_name = "plotfiles";
    experimental::filesystem::create_directory("plotfiles");
    for (int i = 0; i < supernovae.size(); ++i){
        for (int j = 0; j < req_sn_bands.size(); ++j){
            filename = "plotfiles/" + supernovae[i].sn_name + "_" + req_sn_bands[j] + "_band.dat";
            cout << "Opening file: " << filename << " for write" << endl;
            ofstream plotfile(filename);
	    // Use of the '#' char will be ignored by gnuplot
            plotfile << '#' << supernovae[i].sn_name << " " << req_sn_bands[j] << " band light curve" << endl;
            if (s == "lums") plotfile << "#columns: time[days] (relative to peak), Luminosity[erg/s]" << "\n" << endl;
            else plotfile << "#columns: time[days] (relative to peak), Absolute Magnitude[]" << "\n" << endl;
            for (int k = 0; k < supernovae[i].time[j].size(); ++k){
      	        plotfile << supernovae[i].time[j][k] << "  " << supernovae[i].lum[j][k] << endl;
            }      
            cout << "Closing file: " << filename << endl;
            plotfile.close();
        }
    }
  
    return 0;
}    
