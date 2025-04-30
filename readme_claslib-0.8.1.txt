# --------------------------------------
# CLASLIB (CLAS test library)-P/PS
# --------------------------------------

# Overview
CLAS test library is a reference implementation of CLAS, and it contains
conversion utilities or a post-processing PPP-RTK positioning utility depending
on its version. Please note it doesn't assure anything in terms of its 
effectiveness or reliability.

CLAS is the Centimeter-Level Augmentation Service of QZSS.
The service is provided as CLAS messages (Ref. Note 1), the format is
defined in RTCM STANDARD 10403.2(MT4073) and called 'Compact SSR'(Ref. Note 2).

The library is derived from RTKLIB (version 2.4.2p13, Ref. Note 3) and
GSILIB (version.1.0.3, Ref. Note 4).

# Contents of the Library
## CLASLIB Version 0.4.0 and later
The version 0.4.0 and later contains SSR2OSR, which is a Compact SSR to
OSR converter. It also includes a manual and a sample data set. The utilities
DUMPCSSR and QZS2RTCM included in version 0.3.1 are no longer supported from
this version.

The following files are used to decode Compact SSR.
  - cssr.c: functions to decode Compact SSR messages
  - cssr.h: header file for Compact SSR related functions
  - cssr2osr.c: functions for SSR to OSR conversion
  - grid.c: functions to handle the gridded correction

## CLASLIB Version 0.6.0 and later
The version 0.6.0 and later contains SSR2OBS, which is a Compact SSR to
virtual observations converter for the VRS-RTK positioning mode of RNX2RTKP.

## CLASLIB Version 0.7.0 and later
The version 0.7.0 and later supports Compact SSR SubType12 messages 
introduced from IS-QZSS-L6-002.

## CLASLIB Version 0.4.0-P, 0.4.0-PS and later
The version 0.4.0-P and later contains binary executable of RNX2RTKP, which
is a post-processing PPP-RTK positioning utility using compact SSR.
It also includes a manual and a sample data set. The utilities DUMPCSSR and
QZS2RTCM included in version 0.3.1 are no longer supported from this version.

The version 0.4.0-PS and later contains binary executable of RNX2RTKP and
its source code as well. Anything else is the same as version 0.4.0-P.

## CLASLIB Version 0.6.0-P, 0.6.0-PS and later
The version 0.6.0-P and later supports a VRS-RTK positioning mode in RNX2RTKP

It is required to submit an application form to get its source code and/or its
executable. For information on how to submit the application, please refer
to the following WEB page or contact the following e-mail address:

 URL : https://sys.qzss.go.jp/dod/en/downloads/clas.html
 e-mail : IS.PS-QZSS-L6@rm.MitsubishiElectric.co.jp


# Usage:
## SSR2OSR (for version 0.4.0 and later)
To make the executable, move to the folder 'util/ssr2osr' and use 'make' as
follows in a CYGWIN shell.
    > make

To make corrections represented in observation space from sample data,
    > make test1

## SSR2OSR (for version 0.6.1 and later)
To parse l6 message data into each compact ssr sub type message,
    > make test1d
    
In the above command, the option '-dump' is specified.
When this option is specified, only parse process is performed. No other process is performed.

## SSR2OSR (for version 0.7.0 and later)

To make corrections represented in observation space from sample data corresponding to IS-QZSS-L6-002 (including SubType12),
    > make test1_ST12

Please refer to the CLASLIB manual for more information.

## SSR2OBS (for version 0.6.0 and later)
To make the executable, move to the folder 'util/ssr2obs' and use 'make' as
follows in a CYGWIN shell.
    > make

To make virtual observations in RINEX3 format from sample data,
    > make test1r

To make virtual observations in RTCM3 MSM format from sample data,
    > make test1m

## SSR2OBS (for version 0.6.1 and later)
To parse l6 message data into each compact ssr sub type message,
    > make test1d

In the above command, the option '-dump' is specified.
When this option is specified, only parse process is performed. No other process is performed.

Please refer to the CLASLIB manual for more information.

## RNX2RTKP (for version 0.4.0-P/0.4.0-PS and later)
To make the executable, move to the folder 'util/rnx2rtkp' and use 'make' as
follows in a CYGWIN shell. (for version 0.4.0-PS and later only)
    > make

To perform PPP-RTK positioning calculation using sample data,
move to the folder 'util/rnx2rtkp' and do the following in a CYGWIN shell.
    > make test_L6

The results of PPP-RTK positioning can be obtained in the NMEA183-GGA format.

## RNX2RTKP (for version 0.6.0-P/0.6.0-PS and later)
To make the executable, move to the folder 'util/rnx2rtkp' and use 'make' as
follows in a CYGWIN shell. (for version 0.6.0-PS and later only)
    > make

To perform VRS-RTK positioning calculation using sample data,
move to the folder 'util/rnx2rtkp' and do the following in a CYGWIN shell.
    > make test_VRS

## RNX2RTKP (for version 0.6.1-P/0.6.1-PS and later)
To perform positioing caluclation using BINEX format observation data, move to the folder 
'util/rnx2rtkp' and do the following in a CYGWIN shell.
    > test_L6_bnx

For analysis a period around start/end of week epoch,
time consistency check option '-l6w' should be specified.
As the sample, move to the folder 'util/rnx2rtkp' and refer the following in a CYGWIN shell.
    > make test_L6_week

## RNX2RTKP (for version 0.7.0-P/0.7.0-PS and later)
To perform positioing caluclation  from sample data corresponding to IS-QZSS-L6-002 (including SubType12), move to the folder 
'util/rnx2rtkp' and do the following in a CYGWIN shell.
    > test_ST12   
    > test_ST12_VRS

To perform positioing caluclation using BINEX format observation and navigation data, move to the folder 
'util/rnx2rtkp' and do the following in a CYGWIN shell.
    > test_NONAV


## SSR2OBS (for version 0.7.0 and later) 
To make virtual observations in RINEX3 format from sample data corresponding IS-QZSS-L6-002 (including SubType12),
    > make test1r_ST12
    
Please refer to the CLASLIB manual for more information.

# System Requirement:
The executable binaries included in this package require Microsoft Windows 7
64bit system (Ref. Note 5). To compile the CLAS test library, MinGW(Win32) or
CYGWIN is recommended. Other operation systems and compilers can be used with
this library but NOT guaranteed.

# License:
## CLASLIB Version 0.6.0 and earlier
CLASLIB version 0.6.0 and earlier is distributed under the following
BSD 2-clause license (http://opensource.org/licenses/BSD-2-Clause) and
additional two exclusive clauses. Users are permitted to develop, produce
or sell their own non-commercial or commercial products utilizing,
linking or including CLAS Test Library as long as they comply with the license.

 - Copyright (c) 2007-, T. Takasu, All rights reserved.
 - Copyright (c) 2014-, by Geospatial Information Authority of Japan,
All rights reserved.
 - Copyright (c) 2017-, Mitsubishi Electric Corp., All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

- The software package includes some companion executive binaries or shared
  libraries necessary to execute APs on Windows. These licenses succeed to the
  original ones of these software.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Update History
* 2017/09/12   0.2.0 - The first public version
* 2017/12/01   0.3.1 - Fixed several bugs
* 2018/04/10   0.4.0 - Removed conversion utilities: DUMPCSSR and QZS2RTCM.
					   Fixed an issue in compact SSR file reading.
* 2018/04/10   0.4.0-P/0.4.0-PS - Added PPP-RTK post-processing utility (RNX2RTKP).
* 2018/04/23   0.4.1-P/0.4.1-PS - Updated license agreement of 0.4.0-P/0.4.0-PS.
* 2018/05/17   0.4.2 - Updated reference RTKLIB version in readme_claslib.txt
* 2018/08/31   0.5.0 - Supported GPS L5 signal and Galileo
                       Added release notes "relnotes_0.5.0.txt" 
* 2018/09/04   0.5.1 - Fixed an issue in compact SSR invalid value handling
                       Added release notes "relnotes_0.5.1.txt" 
* 2018/09/12   0.5.2 - Added release notes "relnotes_0.5.2.txt" 
* 2019/02/27   0.6.0 - Added a new conversion utility: SSR2OBS 
					 - Added a VRS-RTK positioning mode: RNX2RTKP
 	   				 - Added release notes "relnotes_0.6.0.txt" 
* 2019/09/05   0.6.1 - Supported Binex formart : RNX2RTKP 
					 - Added a dump option: SSR2OSR, SSR2OBS
 	   				 - Added release notes "relnotes_0.6.1.txt" 
* 2019/12/25   0.7.0 - Supported Compact SSR SubType 12 Message
					   : RNX2RTKP, SSR2OSR, SSR2OBS
                     - Supported BINEX upgraded Galileo decoded ephemeris (0x01-14)
                                           : RNX2RTKP, SSR2OSR, SSR2OBS
 	   				 - Added release notes "relnotes_0.7.0.txt
* 2020/02/13   0.7.1 - Refactoring source codes regarding inputs of each cssr subtype message (cssr.c)
 	   				 - Added release notes "relnotes_0.7.1.txt"
* 2020/11/17   0.7.2 - Refactoring source codes and improved fucntions ragarding ocean loading effect
 	   				 - Added release notes "relnotes_0.7.2.txt"
* 2022/08/25   0.7.3 - Improved parameters for PAR and ionospheric residual estimation to optimize for ionospheric disturbances.
 	   				 - Fixed some VRS processes to support PAR.
 	   				 - Added release notes "relnotes_0.7.3.txt"
* 2024/03/11   0.7.3a - There are no changes to SSR2OSR, SSR2OBS and RNX2RTKP.
                        This revision is due to changes of LICENSE.
* 2024/10/03   0.7.4 - Improved several processing.
 	   				 - Added release notes "relnotes_0.7.4.txt"
* 2024/12/06   0.7.5 - There are no changes to SSR2OSR and SSR2OBS.
 	   				   This revision is due to changes in the recommended parameters of the RNX2RTKP configuration files.
 	   				 - Added release notes "relnotes_0.7.5.txt"
* 2025/01/31   0.8.0 - Supported two channels of L6 message input defined in IS-QZSS-L6-006 in ppp-rtk mode
					   : RNX2RTKP, SSR2OSR, SSR2OBS
 	   				 - Added release notes "relnotes_0.8.0.txt"
* 2025/04/30   0.8.1 - Supported two channels of L6 message input defined in IS-QZSS-L6-006 in vrs-rtk mode.
					   : RNX2RTKP, SSR2OBS
					 - Supported settings that prioritizes the L1C signal over the L1C/A signal for QZS in misc-rnxopt1 option in ppp-rtk mode.
					   : RNX2RTKP, SSR2OSR
 	   				 - Added release notes "relnotes_0.8.1.txt"
# Notes
* Note 1: CLAS message archives is avaiable at
  https://sys.qzss.go.jp/dod/archives/clas.html
* Note 2: Compact SSR is defined in the interface specification: IS-QZSS-L6.
  IS-QZSS-L6 is available at
  https://qzss.go.jp/en/technical/ps-is-qzss/ps-is-qzss.html
* Note 3: The RTKLIB is available at GitHub:
  https://github.com/tomojitakasu/RTKLIB.git
* Note 4: The GSILIB is available at GIS web page:
  http://datahouse1.gsi.go.jp/gsilib/gsilib.html
* Note 5: Microsoft(R), Windows(R) are either registered trademarks or
  trademarks of Microsoft Corporation in the United States and/or other countries.
