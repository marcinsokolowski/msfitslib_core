#ifndef _FFT_H__
#define _FFT_H__

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>

// std :
#include <vector>
#include <string>
#include <complex>

using namespace std;

#define N_FINE_CH_PER_COARSE 128
#define N_FINE_CH_PER_BAND   3072 // 24 coarse channels of 128 fine channels 
#define MWA_CLOCK_HZ 655360000

class CBgFits;

class CSpectrometer
{
public :
   CSpectrometer();
   ~CSpectrometer();
   

   static string gPfbCoeffFile;
   static vector<double> gPfbCoeffs;
   static int m_DumpChannel;
   static int m_PolsInFile;
   static int m_Pol;
   static int m_nBits;
   static int m_MaxBytesToProcess;
   static string m_szVoltageDumpFile;
   
   // eda parameters
   static double m_EDA_ElectricalLenM;

   static int doFFT( unsigned char* data_fft, int in_count, double* spectrum, double* spectrum_re, double* spectrum_im, int& out_count, double norm );
   static int doFFT( double* in, int in_count, double* spectrum, double* spectrum_re, double* spectrum_im, int& out_count, double norm );
   static int doFFT( std::complex<float>* in, int in_count, double* spectrum, std::complex<float>* spectrum_reim, int& out_count, double norm );
   static int fileFFT( const char* binfile, double* acc_spec, const char* out_bin_file=NULL, int out_coarse_channel=-1, int skip_extra=0, const char* out_power_file=NULL, const char* out_power_fits=NULL, 
                       int n_out_channels = N_FINE_CH_PER_BAND, time_t file_ux_start=0, long int infile_size_bytes=-1, const char* out_float_file=NULL );
   static int filePFB( const char* binfile, double* acc_spec, const char* out_bin_file, int out_coarse_channel, int skip_extra, int n_taps=12 );


   // binary file reading :
   static int dumpSignatecBinFile( const char* binfile, int dump_idx=0, int n_samples=8192 );      
   static int SkipNSec_and_SaveMSec( const char* binfile, const char* out_bin_file, double skip_n_sec=0, double save_n_sec=2, int skip_extra=0 );
   static int SkipNSpectraAndSaveMSpectra( const char* binfile, const char* out_bin_file, int skip_n_spectra, int save_n_spectra, int spectrum_size=32768, int bytes_per_channel=2 );
   
   static int IntegrateFFT_Power( const char* binfile, vector<double>& avg_spectrum, int spectrum_size=32768, int bytes_per_channel=2 );
   static int CorrelateBinary( const char* eda_file, const char* bighorns_file, 
                               vector<double>& avg_power, vector<double>& avg_eda_power, 
                               vector<double>& cross_power_re, vector<double>& cross_power_im,
                               int spectrum_size=32768, int bytes_per_channel=2,
                               CBgFits* pCrossPowerFullTimeRes=NULL );
   
   static int CorrelateBinaryFloat( const char* eda_file, const char* bighorns_file, 
                               vector<double>& avg_power, vector<double>& avg_eda_power, 
                               vector<double>& cross_power_re, vector<double>& cross_power_im,
                               int spectrum_size=32768, 
                               CBgFits* pCrossPowerFullTimeRes=NULL );
   
   // 
   static int m_DebugNSpectra;
};

#endif
