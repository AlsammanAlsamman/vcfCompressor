//############################################# About Author #########################################
// Created by: Alsamman M. Alsamman                                                                  #
// Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
// License: MIT License - https://opensource.org/licenses/MIT                                        #
// Disclaimer: The script comes with no warranty, use at your own risk                               #
// This script is not intended for commercial use                                                    #
//####################################################################################################


// define these if not defined
#ifndef MAX_WORD_SIZE
#define MAX_WORD_SIZE 50
#endif

#define bool int
#define true 1
#define false 0
#define NA -1 // missing data
#define empty 0 // empty data
// for development purposes
#define BeSilent 1

// printing to the log file
#define LOG(...) do{FILE *LOGFILE = fopen("log.txt", "a"); fprintf(LOGFILE, __VA_ARGS__); fclose(LOGFILE); if(BeSilent) (__VA_ARGS__);} while (0)
// printing to the error file
#define ERR(...) do{FILE *ERRFILE = fopen("error.txt", "a"); fprintf(ERRFILE, __VA_ARGS__); fclose(ERRFILE); if(BeSilent) printf(__VA_ARGS__);} while (0)
// EXIT ERROR
#define EXITERR(...) ERR(__VA_ARGS__); exit(1);