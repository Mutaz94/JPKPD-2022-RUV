Sun Oct 24 02:53:26 CDT 2021
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
$DATA ../../../../data/SD4/S1/dat98.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       24 OCT 2021
Days until program expires : 175
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
 NO. OF DATA RECS IN DATA SET:      500
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

 TOT. NO. OF OBS RECS:      400
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m98.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1681.59732232564        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6539E+02 -1.7029E+01 -3.5962E+01  5.3628E+01  6.1370E+01  6.1948E+01 -4.3231E+00  2.6103E+00  2.0077E+01  4.2337E+00
             4.9130E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1684.90234307993        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      132
 NPARAMETR:  9.6833E-01  1.0790E+00  1.0442E+00  1.0070E+00  9.4555E-01  8.8350E-01  1.0175E+00  9.9592E-01  9.6226E-01  9.8174E-01
             8.9520E-01
 PARAMETER:  6.7815E-02  1.7607E-01  1.4324E-01  1.0701E-01  4.4013E-02 -2.3862E-02  1.1737E-01  9.5912E-02  6.1532E-02  8.1575E-02
            -1.0711E-02
 GRADIENT:  -1.2231E+01  6.1587E+01  3.3976E+01  4.7323E+01 -8.0929E+01 -2.5951E+01 -8.2373E+00 -5.1651E+00  3.8549E+00  2.8932E+00
             5.3872E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1687.72659945791        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      308
 NPARAMETR:  9.7230E-01  1.0851E+00  1.2308E+00  1.0017E+00  1.0356E+00  9.2748E-01  1.2632E+00  1.4592E+00  9.1210E-01  9.2075E-01
             8.6086E-01
 PARAMETER:  7.1910E-02  1.8163E-01  3.0763E-01  1.0171E-01  1.3501E-01  2.4717E-02  3.3368E-01  4.7789E-01  7.9980E-03  1.7436E-02
            -4.9820E-02
 GRADIENT:   2.4478E+00  4.1573E+01  1.5229E+01  3.0948E+01 -2.6531E+01 -5.0071E+00  1.2435E+01 -7.7205E-01  5.5109E+00 -5.4846E+00
            -8.4380E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1690.35514129246        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      488
 NPARAMETR:  9.7186E-01  1.2214E+00  8.8491E-01  8.7603E-01  1.0034E+00  9.4101E-01  1.0580E+00  1.0849E+00  9.5411E-01  9.3011E-01
             8.7576E-01
 PARAMETER:  7.1459E-02  2.9998E-01 -2.2275E-02 -3.2357E-02  1.0340E-01  3.9194E-02  1.5642E-01  1.8147E-01  5.3019E-02  2.7551E-02
            -3.2663E-02
 GRADIENT:  -4.2546E+00  1.7844E+00 -3.4102E+00  8.1063E+00  3.4553E+00 -5.2601E-02  1.1411E+00  6.4754E-01  7.1970E-01  2.5650E-01
            -4.5347E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1690.57875308630        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      663
 NPARAMETR:  9.7565E-01  1.4708E+00  7.1552E-01  7.1096E-01  1.0500E+00  9.4216E-01  9.2184E-01  1.0198E+00  1.0819E+00  9.4930E-01
             8.7751E-01
 PARAMETER:  7.5346E-02  4.8583E-01 -2.3475E-01 -2.4114E-01  1.4882E-01  4.0424E-02  1.8613E-02  1.1963E-01  1.7876E-01  4.7975E-02
            -3.0667E-02
 GRADIENT:   2.7276E+00  6.0788E+00  3.4447E+00  1.1906E+00 -5.4040E+00  8.2312E-02 -5.9288E-01 -7.2398E-01 -7.6715E-01 -3.0920E-01
            -1.0997E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1690.62613148072        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      843
 NPARAMETR:  9.7546E-01  1.5265E+00  6.6231E-01  6.6994E-01  1.0612E+00  9.4236E-01  8.9971E-01  1.0179E+00  1.1213E+00  9.5374E-01
             8.7732E-01
 PARAMETER:  7.5151E-02  5.2295E-01 -3.1202E-01 -3.0056E-01  1.5941E-01  4.0637E-02 -5.6783E-03  1.1771E-01  2.1451E-01  5.2634E-02
            -3.0879E-02
 GRADIENT:   1.6884E+00 -2.8679E+00  3.8427E-01 -4.5040E-01 -3.2178E-01  8.2698E-02 -1.6094E-02 -9.7270E-02  1.6278E-02  4.7559E-01
             1.8511E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1690.62826794577        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1019
 NPARAMETR:  9.7534E-01  1.5277E+00  6.6083E-01  6.6934E-01  1.0609E+00  9.4240E-01  8.9957E-01  1.0305E+00  1.1212E+00  9.4860E-01
             8.7689E-01
 PARAMETER:  7.5034E-02  5.2374E-01 -3.1425E-01 -3.0146E-01  1.5915E-01  4.0672E-02 -5.8432E-03  1.3002E-01  2.1436E-01  4.7232E-02
            -3.1369E-02
 GRADIENT:   1.3731E+00 -2.7075E+00  1.3242E-01  1.9292E-02  4.4061E-01  9.2137E-02  2.0728E-02 -9.2607E-03  7.5962E-03 -3.2927E-02
            -1.1427E-02

0ITERATION NO.:   31    OBJECTIVE VALUE:  -1690.62826794577        NO. OF FUNC. EVALS.:  24
 CUMULATIVE NO. OF FUNC. EVALS.:     1043
 NPARAMETR:  9.7534E-01  1.5277E+00  6.6083E-01  6.6934E-01  1.0609E+00  9.4240E-01  8.9957E-01  1.0305E+00  1.1212E+00  9.4860E-01
             8.7689E-01
 PARAMETER:  7.5034E-02  5.2374E-01 -3.1425E-01 -3.0146E-01  1.5915E-01  4.0672E-02 -5.8432E-03  1.3002E-01  2.1436E-01  4.7232E-02
            -3.1369E-02
 GRADIENT:  -6.1979E-01  5.9154E-01  1.3837E-01  3.2533E-01  4.7255E-01 -7.9012E-02  1.2232E-02 -5.1204E-03 -1.1863E-02 -3.4137E-02
            -1.0980E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1043
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.5738E-04 -2.0573E-02 -2.8583E-02  1.7531E-02 -3.6530E-02
 SE:             2.9890E-02  2.3929E-02  1.0891E-02  2.1979E-02  2.1649E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8779E-01  3.8992E-01  8.6780E-03  4.2508E-01  9.1535E-02

 ETASHRINKSD(%)  1.0000E-10  1.9834E+01  6.3514E+01  2.6369E+01  2.7473E+01
 ETASHRINKVR(%)  1.0000E-10  3.5734E+01  8.6688E+01  4.5785E+01  4.7398E+01
 EBVSHRINKSD(%)  3.6286E-01  1.9812E+01  6.7092E+01  2.7997E+01  2.4587E+01
 EBVSHRINKVR(%)  7.2441E-01  3.5699E+01  8.9170E+01  4.8155E+01  4.3130E+01
 RELATIVEINF(%)  9.9101E+01  2.7160E+00  6.9635E-01  1.9858E+00  1.1034E+01
 EPSSHRINKSD(%)  4.6082E+01
 EPSSHRINKVR(%)  7.0929E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1690.6282679457663     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -955.47744138202813     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     4.74
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1690.628       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.75E-01  1.53E+00  6.61E-01  6.69E-01  1.06E+00  9.42E-01  9.00E-01  1.03E+00  1.12E+00  9.49E-01  8.77E-01
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,       29.885
Stop Time:
Sun Oct 24 02:53:33 CDT 2021
