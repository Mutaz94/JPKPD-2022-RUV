Wed Sep 29 21:52:00 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	Two-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa1/A1/dat14.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
$PK
ET1 = EXP(ETA(1)*THETA(6))
ET2 = EXP(ETA(2)*THETA(7))
ET3 = EXP(ETA(3)*THETA(8))
ET4 = EXP(ETA(4)*THETA(9))
ET5 = EXP(ETA(5)*THETA(10))

CL = 5.0 * THETA(1) * ET1
V2 = 35  * THETA(2) * ET2
Q  = 50  * THETA(3) * ET3
V3 = 50  * THETA(4) * ET4
KA = 0.7 * THETA(5) * ET5
SC = V2
$ERROR
CVERR = 0.05
W = THETA(11)*F*CVERR

Y 	= F + W*ERR(1)

$THETA
(0,1) ; CL
(0,1) ; V2
(0,1) ; Q
(0,1) ; V3
(0,1) ; KA
(0,1) ; IIVCL
(0,1) ; IIVV2
(0,1) ; IIVQ
(0,1) ; IIVV3
(0,1) ; IIVKA
(0,1) ; CVPropErr

$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       29 SEP 2021
Days until program expires : 200
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 template control stream
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      600
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      500
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:  11
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 TWO COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN4)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K23)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K32)
   BASIC PK PARAMETER NO.  5: ABSORPTION RATE (KA)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V2, Q, V3 TO K, K23, K32 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         PERIPH.      ON         NO         YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            6           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
1


 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            10000
 NO. OF SIG. FIGURES REQUIRED:            2
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): m14.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1709.76499563253        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2068E+02 -2.3189E+01 -2.2900E+01  3.2040E+01  1.3926E+02  3.2425E+01 -2.0402E+01  6.8117E+00  2.8890E+01 -3.3958E+01
            -7.7256E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1892.58539654351        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0170E+00  1.0333E+00  1.0367E+00  1.0460E+00  9.1766E-01  1.0774E+00  1.0467E+00  8.4234E-01  8.0642E-01  9.6191E-01
             1.6568E+00
 PARAMETER:  1.1689E-01  1.3275E-01  1.3603E-01  1.4501E-01  1.4073E-02  1.7456E-01  1.4561E-01 -7.1567E-02 -1.1514E-01  6.1165E-02
             6.0490E-01
 GRADIENT:   1.9623E+02  4.9043E+01  1.7912E+01  5.7432E+01 -1.9434E+01  4.5545E+01 -9.7250E+00  5.7062E+00 -1.0813E+01 -4.7395E+00
            -2.7231E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1900.52951025115        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      181
 NPARAMETR:  1.0078E+00  8.5735E-01  6.3585E-01  1.1104E+00  6.7516E-01  1.0289E+00  1.3335E+00  2.1289E-01  7.6779E-01  8.3806E-01
             1.6131E+00
 PARAMETER:  1.0778E-01 -5.3914E-02 -3.5279E-01  2.0474E-01 -2.9280E-01  1.2854E-01  3.8780E-01 -1.4470E+00 -1.6424E-01 -7.6663E-02
             5.7819E-01
 GRADIENT:  -7.1399E+00  1.3323E+01 -3.6605E+01  4.6132E+01  3.8667E+01 -9.3941E-01  4.2924E-01  1.2799E+00 -4.7436E+00  1.1505E+01
            -3.0331E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1906.36533353118        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      358
 NPARAMETR:  1.0071E+00  5.8255E-01  6.4353E-01  1.2413E+00  5.7227E-01  1.0197E+00  1.6795E+00  6.0010E-02  7.4908E-01  7.5142E-01
             1.6948E+00
 PARAMETER:  1.0712E-01 -4.4034E-01 -3.4078E-01  3.1616E-01 -4.5815E-01  1.1953E-01  6.1851E-01 -2.7132E+00 -1.8891E-01 -1.8579E-01
             6.2758E-01
 GRADIENT:  -3.7781E+00  9.3102E+00  5.1685E-01  7.9872E+00 -1.1155E+01 -2.6285E+00  6.5655E-01  1.0370E-01 -6.8873E-01  4.6509E+00
             1.5079E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1907.52945697586        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      533
 NPARAMETR:  1.0071E+00  4.1993E-01  7.9722E-01  1.3504E+00  6.1392E-01  1.0215E+00  2.0502E+00  2.9746E-02  7.3705E-01  8.6418E-01
             1.6611E+00
 PARAMETER:  1.0706E-01 -7.6766E-01 -1.2663E-01  4.0042E-01 -3.8790E-01  1.2124E-01  8.1795E-01 -3.4151E+00 -2.0510E-01 -4.5978E-02
             6.0746E-01
 GRADIENT:   6.9793E+00  3.6077E+00  5.5111E+00 -6.1487E+00 -1.0462E+01 -4.0955E-01  1.7934E+00  1.8916E-02 -9.2768E-03  3.5402E+00
            -2.5300E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1908.35151508382        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      711
 NPARAMETR:  9.9414E-01  2.4366E-01  9.5402E-01  1.4895E+00  6.5113E-01  1.0162E+00  2.3979E+00  1.0000E-02  7.4437E-01  9.3720E-01
             1.6833E+00
 PARAMETER:  9.4122E-02 -1.3120E+00  5.2934E-02  4.9841E-01 -3.2904E-01  1.1609E-01  9.7459E-01 -5.0739E+00 -1.9521E-01  3.5145E-02
             6.2078E-01
 GRADIENT:  -1.0947E+01  5.3611E+00  8.4451E-01  3.0358E+01 -5.0101E-01 -1.3902E+00  3.7098E-01  0.0000E+00  4.3010E+00 -1.1854E+00
             4.4517E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1909.50751062479        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      890
 NPARAMETR:  9.9485E-01  7.7309E-02  9.6246E-01  1.5742E+00  6.1813E-01  1.0130E+00  3.2161E+00  1.0000E-02  7.2278E-01  9.5384E-01
             1.6747E+00
 PARAMETER:  9.4840E-02 -2.4599E+00  6.1741E-02  5.5377E-01 -3.8105E-01  1.1293E-01  1.2682E+00 -9.8927E+00 -2.2465E-01  5.2737E-02
             6.1562E-01
 GRADIENT:  -3.2845E-01  1.1712E+00  3.0872E+00  2.1365E+01 -7.2340E+00 -1.4919E+00 -9.0198E-01  0.0000E+00  1.5103E+00  1.0248E+00
             5.5613E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1910.10837841455        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1066
 NPARAMETR:  9.9310E-01  1.0053E-02  9.6173E-01  1.6029E+00  6.0522E-01  1.0175E+00  5.6389E+00  1.0000E-02  7.0825E-01  9.5136E-01
             1.6716E+00
 PARAMETER:  9.3074E-02 -4.4999E+00  6.0980E-02  5.7184E-01 -4.0216E-01  1.1737E-01  1.8297E+00 -1.9100E+01 -2.4495E-01  5.0133E-02
             6.1377E-01
 GRADIENT:   1.3639E-01  2.0665E-01  2.6371E+00  5.4598E+00 -5.6747E+00  7.3832E-01 -8.1791E-02  0.0000E+00 -1.0845E+00  5.6923E-01
            -7.1864E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1910.14266003197        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1257
 NPARAMETR:  9.9317E-01  1.0000E-02  9.6154E-01  1.5979E+00  6.0585E-01  1.0157E+00  6.1155E+00  1.0000E-02  7.1059E-01  9.4999E-01
             1.6725E+00
 PARAMETER:  9.3142E-02 -4.5603E+00  6.0779E-02  5.6868E-01 -4.0113E-01  1.1559E-01  1.9108E+00 -1.9100E+01 -2.4165E-01  4.8700E-02
             6.1432E-01
 GRADIENT:   4.5085E-01  0.0000E+00  2.5812E+00 -7.6230E+00 -3.0950E+00  1.1188E-01 -8.6535E-02  0.0000E+00  1.2863E-01  2.2184E-01
             3.3119E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1910.15669031678        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1441
 NPARAMETR:  9.9304E-01  1.0000E-02  9.5942E-01  1.5995E+00  6.0605E-01  1.0155E+00  6.9637E+00  1.0000E-02  7.1050E-01  9.4862E-01
             1.6724E+00
 PARAMETER:  9.3013E-02 -4.5603E+00  5.8578E-02  5.6969E-01 -4.0079E-01  1.1538E-01  2.0407E+00 -1.9100E+01 -2.4178E-01  4.7257E-02
             6.1426E-01
 GRADIENT:   1.7788E-01  0.0000E+00 -2.2684E-02 -2.8419E+00 -1.2953E-01  1.3642E-02 -9.8716E-02  0.0000E+00  7.7091E-02  3.7488E-04
            -6.0288E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1910.17188779064        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     1606
 NPARAMETR:  9.9267E-01  1.0000E-02  9.6153E-01  1.6051E+00  6.0677E-01  1.0154E+00  1.1347E+01  1.0000E-02  7.0939E-01  9.4872E-01
             1.6731E+00
 PARAMETER:  9.2646E-02 -4.5603E+00  6.0773E-02  5.7321E-01 -3.9960E-01  1.1526E-01  2.5290E+00 -1.9100E+01 -2.4335E-01  4.7356E-02
             6.1468E-01
 GRADIENT:  -7.8905E-01  0.0000E+00 -7.1152E-01  1.0368E+01 -1.1167E+00 -7.8922E-02 -1.0492E-02  0.0000E+00 -1.8816E-01  2.3481E-02
             8.9152E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1606
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2452E-03  1.1203E-03 -7.6028E-05 -1.1776E-02 -1.6440E-02
 SE:             2.9673E-02  2.1495E-03  1.7222E-04  2.8584E-02  2.4289E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6653E-01  6.0222E-01  6.5887E-01  6.8035E-01  4.9849E-01

 ETASHRINKSD(%)  5.9211E-01  9.2799E+01  9.9423E+01  4.2402E+00  1.8630E+01
 ETASHRINKVR(%)  1.1807E+00  9.9481E+01  9.9997E+01  8.3005E+00  3.3790E+01
 EBVSHRINKSD(%)  8.4317E-01  9.3652E+01  9.9378E+01  4.1296E+00  1.7423E+01
 EBVSHRINKVR(%)  1.6792E+00  9.9597E+01  9.9996E+01  8.0887E+00  3.1811E+01
 RELATIVEINF(%)  9.3111E+01  1.6492E-02  2.9893E-04  4.4628E+00  6.1113E+00
 EPSSHRINKSD(%)  2.9468E+01
 EPSSHRINKVR(%)  5.0252E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1910.1718877906387     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -991.23335458596603     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.05
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.56
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1910.172       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.93E-01  1.00E-02  9.62E-01  1.61E+00  6.07E-01  1.02E+00  1.13E+01  1.00E-02  7.09E-01  9.49E-01  1.67E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  9.00E-02
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  3.00E-01
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.00E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.08E+03
 
 TH 2
+        0.00E+00  3.00E+03
 
 TH 3
+       -9.46E+00  0.00E+00  4.69E+02
 
 TH 4
+       -1.40E+01  0.00E+00 -1.19E+02  7.90E+02
 
 TH 5
+        1.74E+01  0.00E+00 -9.44E+02 -7.60E+01  2.30E+03
 
 TH 6
+        2.59E+00  0.00E+00  6.77E-01 -4.81E+00 -1.33E+00  1.87E+02
 
 TH 7
+        6.46E-03  0.00E+00 -2.80E-02 -1.08E-01  9.21E-02  4.83E-04  2.52E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.33E+00  0.00E+00  1.99E+01 -1.88E+01 -5.36E+00  1.47E+00 -3.34E-02  0.00E+00  3.34E+02
 
 TH10
+        1.04E+00  0.00E+00 -1.36E+01  5.12E+00 -7.01E+01  4.10E-01 -5.06E-02  0.00E+00  2.29E-01  1.05E+02
 
 TH11
+       -1.12E+01  0.00E+00 -1.33E+01 -1.46E+01  1.29E+01  1.15E+00 -1.55E-02  0.00E+00  1.07E+01  2.28E+01  1.58E+02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       31.686
Stop Time:
Wed Sep 29 21:52:34 CDT 2021
