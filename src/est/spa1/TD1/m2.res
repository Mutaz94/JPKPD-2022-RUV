Thu Sep 30 01:05:21 CDT 2021
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
$DATA ../../../../data/spa1/TD1/dat2.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m2.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2144.15578103754        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4845E+02 -2.0922E+01 -1.7383E+01  8.4693E+00 -3.0529E+01  3.6830E+01  1.1027E+01  1.9447E+01  4.4873E+01  1.5049E+01
             8.9509E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2155.55539191294        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.7218E-01  9.9802E-01  1.2195E+00  1.0778E+00  1.1690E+00  9.3433E-01  9.2256E-01  7.8792E-01  7.5092E-01  1.0385E+00
             9.9329E-01
 PARAMETER:  7.1781E-02  9.8019E-02  2.9841E-01  1.7491E-01  2.5616E-01  3.2078E-02  1.9393E-02 -1.3836E-01 -1.8645E-01  1.3780E-01
             9.3269E-02
 GRADIENT:   3.9194E+02  4.0855E+01 -3.0302E+01  1.9198E+02  6.2648E+01  1.2942E+01 -8.9071E+00  3.7252E+00 -9.6963E+00 -1.2970E+01
            -3.8681E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2157.38076741513        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      203
 NPARAMETR:  9.6513E-01  1.0101E+00  1.3143E+00  1.0641E+00  1.2302E+00  9.4197E-01  9.6820E-01  5.8856E-01  7.9544E-01  1.1793E+00
             9.7025E-01
 PARAMETER:  6.4509E-02  1.1004E-01  3.7332E-01  1.6211E-01  3.0716E-01  4.0218E-02  6.7683E-02 -4.3007E-01 -1.2887E-01  2.6491E-01
             6.9798E-02
 GRADIENT:  -3.3927E+01  4.0672E+00 -9.7456E+00  2.7128E+01  3.5242E+01 -2.3588E+01 -2.1042E+00 -2.5649E+00 -5.4711E+00 -3.2949E+00
            -2.5710E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2159.68937471265        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      379
 NPARAMETR:  9.7798E-01  1.0129E+00  1.3335E+00  1.0546E+00  1.1980E+00  9.9549E-01  9.5922E-01  7.5446E-01  8.1431E-01  1.1553E+00
             9.9948E-01
 PARAMETER:  7.7733E-02  1.1279E-01  3.8780E-01  1.5319E-01  2.8066E-01  9.5475E-02  5.8366E-02 -1.8175E-01 -1.0542E-01  2.4434E-01
             9.9476E-02
 GRADIENT:  -1.4965E+00  2.5072E-01 -8.6915E-01  2.2886E+00  6.9018E-01  1.1383E-01 -1.7835E-01 -5.9906E-04 -1.3222E+00  2.8460E-01
             2.0894E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2159.83778202404        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      561
 NPARAMETR:  9.7947E-01  1.1819E+00  1.1818E+00  9.4757E-01  1.2094E+00  9.9344E-01  8.7116E-01  6.3176E-01  8.8602E-01  1.1433E+00
             9.9811E-01
 PARAMETER:  7.9260E-02  2.6715E-01  2.6700E-01  4.6146E-02  2.9012E-01  9.3417E-02 -3.7930E-02 -3.5924E-01 -2.1015E-02  2.3391E-01
             9.8107E-02
 GRADIENT:  -1.5556E+00  9.4807E+00  3.1593E+00  8.5857E+00 -6.2894E+00 -1.2678E+00  5.0659E-01 -1.6609E-01 -1.4490E-01 -2.7768E-01
             8.6353E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2159.98162686658        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      738
 NPARAMETR:  9.8217E-01  1.3234E+00  1.0659E+00  8.5661E-01  1.2390E+00  9.9969E-01  7.9227E-01  5.3560E-01  9.5689E-01  1.1513E+00
             9.9485E-01
 PARAMETER:  8.2006E-02  3.8024E-01  1.6383E-01 -5.4771E-02  3.1427E-01  9.9688E-02 -1.3286E-01 -5.2436E-01  5.5928E-02  2.4091E-01
             9.4833E-02
 GRADIENT:   2.2186E+00  1.2695E+01  1.2733E+00  1.2642E+01 -2.8013E+00  8.0111E-01 -1.0635E+00 -5.8200E-02 -1.0576E+00 -3.0379E-01
            -1.8121E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2160.21075922506        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      915
 NPARAMETR:  9.8351E-01  1.5324E+00  8.8576E-01  7.1750E-01  1.2806E+00  1.0007E+00  7.2926E-01  3.5591E-01  1.0735E+00  1.1580E+00
             9.9672E-01
 PARAMETER:  8.3374E-02  5.2683E-01 -2.1305E-02 -2.3199E-01  3.4729E-01  1.0066E-01 -2.1573E-01 -9.3307E-01  1.7093E-01  2.4673E-01
             9.6710E-02
 GRADIENT:   2.2212E+00  1.2776E+01  1.0585E-01  9.4504E+00 -3.6648E+00  6.9623E-01 -1.1155E+00  1.7758E-01 -1.6518E+00  1.4591E-01
             4.2253E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2160.32382899337        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1092
 NPARAMETR:  9.8285E-01  1.6710E+00  7.9837E-01  6.2592E-01  1.3340E+00  9.9906E-01  6.9183E-01  2.7080E-01  1.1906E+00  1.1829E+00
             9.9697E-01
 PARAMETER:  8.2705E-02  6.1342E-01 -1.2518E-01 -3.6853E-01  3.8822E-01  9.9058E-02 -2.6842E-01 -1.2064E+00  2.7448E-01  2.6794E-01
             9.6967E-02
 GRADIENT:  -4.7830E-01  1.1858E+01  1.6628E-02  8.1075E+00 -2.6787E+00 -1.5071E-01 -6.0998E-01  1.2973E-01 -7.0547E-01 -3.2085E-02
             4.9313E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2160.43208939554        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1273
 NPARAMETR:  9.8306E-01  1.6744E+00  7.9271E-01  6.1764E-01  1.3399E+00  9.9944E-01  6.9182E-01  9.4771E-02  1.2072E+00  1.1889E+00
             9.9625E-01
 PARAMETER:  8.2912E-02  6.1543E-01 -1.3230E-01 -3.8184E-01  3.9257E-01  9.9442E-02 -2.6843E-01 -2.2563E+00  2.8830E-01  2.7301E-01
             9.6242E-02
 GRADIENT:   7.3442E-02 -6.6290E-01  7.3346E-01  1.5510E+00 -8.1798E-01 -8.4610E-04  1.6770E-01  1.1033E-02  7.0627E-02  2.1549E-01
             2.6927E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2160.44615857151        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1460
 NPARAMETR:  9.8427E-01  1.7187E+00  7.6082E-01  5.8444E-01  1.3591E+00  9.9987E-01  6.8066E-01  1.0000E-02  1.2501E+00  1.1962E+00
             9.9622E-01
 PARAMETER:  8.4145E-02  6.4159E-01 -1.7335E-01 -4.3710E-01  4.0682E-01  9.9874E-02 -2.8470E-01 -5.0904E+00  3.2323E-01  2.7913E-01
             9.6212E-02
 GRADIENT:   2.4796E+00 -9.8594E+00  2.2702E-01 -2.6966E+00  9.0015E-02  1.2858E-01  1.2138E-01  0.0000E+00 -1.4162E-01  8.8565E-03
             7.3555E-02

0ITERATION NO.:   47    OBJECTIVE VALUE:  -2160.44900195399        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1517
 NPARAMETR:  9.8329E-01  1.7214E+00  7.6049E-01  5.8625E-01  1.3590E+00  9.9973E-01  6.8037E-01  1.0000E-02  1.2507E+00  1.1961E+00
             9.9618E-01
 PARAMETER:  8.3145E-02  6.4316E-01 -1.7379E-01 -4.3401E-01  4.0678E-01  9.9729E-02 -2.8512E-01 -5.0904E+00  3.2367E-01  2.7908E-01
             9.6173E-02
 GRADIENT:   2.1428E-01 -2.4021E+00  4.1886E-02  1.4331E+00  1.1059E-01  5.2565E-02  4.5534E-02  0.0000E+00  4.2100E-02 -1.8450E-02
            -5.0829E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1517
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0491E-03 -3.1223E-02 -2.7133E-04  2.3231E-02 -3.7902E-02
 SE:             2.9869E-02  2.1613E-02  1.0198E-04  2.2855E-02  2.3434E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7198E-01  1.4856E-01  7.7977E-03  3.0942E-01  1.0579E-01

 ETASHRINKSD(%)  1.0000E-10  2.7594E+01  9.9658E+01  2.3432E+01  2.1493E+01
 ETASHRINKVR(%)  1.0000E-10  4.7573E+01  9.9999E+01  4.1374E+01  3.8367E+01
 EBVSHRINKSD(%)  3.3088E-01  2.5740E+01  9.9704E+01  2.6381E+01  1.8362E+01
 EBVSHRINKVR(%)  6.6067E-01  4.4854E+01  9.9999E+01  4.5803E+01  3.3352E+01
 RELATIVEINF(%)  9.9187E+01  2.5579E+00  1.2276E-04  2.5533E+00  1.7846E+01
 EPSSHRINKSD(%)  3.2425E+01
 EPSSHRINKVR(%)  5.4336E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2160.4490019539930     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1241.5104687493204     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.79
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.31
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2160.449       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.83E-01  1.72E+00  7.60E-01  5.86E-01  1.36E+00  1.00E+00  6.80E-01  1.00E-02  1.25E+00  1.20E+00  9.96E-01
 


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
+        1.14E+03
 
 TH 2
+       -6.10E+00  4.71E+02
 
 TH 3
+        4.39E+00  1.01E+02  2.01E+02
 
 TH 4
+       -6.01E+00  5.12E+02 -1.75E+02  1.13E+03
 
 TH 5
+        1.71E+00 -1.14E+02 -1.52E+02  1.51E+02  2.81E+02
 
 TH 6
+        2.01E+00 -1.29E+00  1.08E+00 -2.71E+00 -5.06E-01  1.97E+02
 
 TH 7
+        1.26E+00 -3.66E+00  9.44E+00 -1.70E+01 -1.58E+01 -4.37E-01  1.38E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.05E-01 -1.94E+01 -2.53E+01  5.46E+01  4.89E+00 -3.18E-01  3.63E+01  0.00E+00  4.64E+01
 
 TH10
+        1.31E+00 -1.02E+01 -2.26E+01  6.07E+00 -4.04E+01  9.22E-02  5.96E+00  0.00E+00  4.75E+00  6.69E+01
 
 TH11
+       -8.19E+00 -2.47E+01 -2.90E+01  2.11E+00  2.74E+00  1.81E+00  9.97E+00  0.00E+00  5.68E+00  1.46E+01  4.25E+02
 
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
 #CPUT: Total CPU Time in Seconds,       30.162
Stop Time:
Thu Sep 30 01:05:53 CDT 2021
