Sat Oct 23 20:35:38 CDT 2021
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
$DATA ../../../../data/SD2/D2/dat17.csv ignore=@
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
Current Date:       23 OCT 2021
Days until program expires : 176
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
 NO. OF DATA RECS IN DATA SET:      800
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

 TOT. NO. OF OBS RECS:      700
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
 RAW OUTPUT FILE (FILE): m17.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2768.36182950430        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.8662E+02 -2.8122E+01  3.4236E+01 -1.1799E+02  3.1929E+01 -8.0360E+01 -7.6383E+01 -1.7066E+01 -1.2387E+02 -5.5654E+01
            -5.8938E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2833.36893919649        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      130
 NPARAMETR:  1.1467E+00  9.4285E-01  8.8219E-01  1.0454E+00  9.0528E-01  1.2993E+00  1.4051E+00  1.0650E+00  1.5514E+00  1.2855E+00
             1.0172E+00
 PARAMETER:  2.3686E-01  4.1150E-02 -2.5345E-02  1.4442E-01  4.9374E-04  3.6179E-01  4.4009E-01  1.6302E-01  5.3913E-01  3.5114E-01
             1.1703E-01
 GRADIENT:   9.7730E+01 -3.0149E+01 -2.7089E+01 -6.6942E+01 -1.9059E+01 -3.1719E+01  5.8241E+00  6.9496E+00  4.0914E+01  1.3097E+01
            -9.7750E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2843.41710250524        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      313
 NPARAMETR:  1.1576E+00  1.0459E+00  1.0517E+00  1.0600E+00  1.1566E+00  1.4159E+00  1.8630E+00  1.3166E+00  1.2510E+00  1.2983E+00
             1.0621E+00
 PARAMETER:  2.4632E-01  1.4487E-01  1.5040E-01  1.5822E-01  2.4546E-01  4.4774E-01  7.2217E-01  3.7502E-01  3.2391E-01  3.6105E-01
             1.6029E-01
 GRADIENT:   9.4691E+01 -3.8407E+01 -3.7080E+01  3.6321E+00  4.9379E+01  7.2541E+00  7.0583E+00  6.0961E+00  2.6665E+01  1.9436E+01
             3.1836E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2852.81409873055        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      489
 NPARAMETR:  1.0741E+00  1.0660E+00  1.0993E+00  1.1003E+00  1.0764E+00  1.3828E+00  1.8477E+00  1.2569E+00  1.1609E+00  1.1307E+00
             1.0378E+00
 PARAMETER:  1.7147E-01  1.6394E-01  1.9465E-01  1.9558E-01  1.7365E-01  4.2409E-01  7.1393E-01  3.2868E-01  2.4923E-01  2.2288E-01
             1.3713E-01
 GRADIENT:   7.9666E+00  7.1444E+00 -9.8458E-01  1.5546E+01 -4.1888E+00  5.8479E+00  2.1475E+00 -5.9644E-01  1.0142E+01  2.1537E+00
            -4.0398E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2853.33749670609        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      667             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0746E+00  1.0765E+00  1.1211E+00  1.0760E+00  1.1056E+00  1.4211E+00  1.8688E+00  1.3043E+00  1.0817E+00  1.1295E+00
             1.0434E+00
 PARAMETER:  1.7194E-01  1.7370E-01  2.1431E-01  1.7323E-01  2.0043E-01  4.5145E-01  7.2527E-01  3.6563E-01  1.7855E-01  2.2180E-01
             1.4246E-01
 GRADIENT:   7.0805E+02  8.7871E+01  1.0436E+01  1.4792E+02  3.8021E+01  4.7988E+02  1.4896E+02  2.7016E+00  1.4952E+01  6.0371E+00
             1.9549E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2853.40714323825        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      806             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0702E+00  1.0789E+00  1.1169E+00  1.0758E+00  1.1062E+00  1.3967E+00  1.8725E+00  1.3036E+00  1.0867E+00  1.1294E+00
             1.0433E+00
 PARAMETER:  1.6786E-01  1.7597E-01  2.1060E-01  1.7305E-01  2.0096E-01  4.3409E-01  7.2725E-01  3.6511E-01  1.8313E-01  2.2169E-01
             1.4241E-01
 GRADIENT:   6.8669E+02  8.9818E+01  8.8434E+00  1.4919E+02  3.9543E+01  4.5721E+02  1.5010E+02  2.9474E+00  1.6291E+01  6.0160E+00
             1.9295E+00

0ITERATION NO.:   29    OBJECTIVE VALUE:  -2853.41130239426        NO. OF FUNC. EVALS.: 132
 CUMULATIVE NO. OF FUNC. EVALS.:      938
 NPARAMETR:  1.0705E+00  1.0783E+00  1.1188E+00  1.0747E+00  1.1058E+00  1.3960E+00  1.8679E+00  1.3009E+00  1.0821E+00  1.1295E+00
             1.0433E+00
 PARAMETER:  1.6813E-01  1.7538E-01  2.1225E-01  1.7204E-01  2.0059E-01  4.3360E-01  7.2480E-01  3.6307E-01  1.7893E-01  2.2182E-01
             1.4244E-01
 GRADIENT:  -3.1006E-01  3.7620E-02  3.7237E-01 -2.3350E-01 -6.3623E-01  2.0002E-01  2.8940E-02 -8.2133E-02 -1.2643E-01 -7.0165E-02
            -2.2625E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      938
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.1625E-04 -6.0320E-03 -4.0542E-02  4.6383E-03 -3.1601E-02
 SE:             2.9921E-02  2.4246E-02  1.7348E-02  2.4958E-02  2.1916E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8357E-01  8.0353E-01  1.9440E-02  8.5257E-01  1.4933E-01

 ETASHRINKSD(%)  1.0000E-10  1.8774E+01  4.1882E+01  1.6388E+01  2.6578E+01
 ETASHRINKVR(%)  1.0000E-10  3.4024E+01  6.6223E+01  3.0091E+01  4.6092E+01
 EBVSHRINKSD(%)  1.5733E-01  1.7077E+01  4.5225E+01  1.9317E+01  2.5197E+01
 EBVSHRINKVR(%)  3.1440E-01  3.1237E+01  6.9998E+01  3.4903E+01  4.4045E+01
 RELATIVEINF(%)  9.9682E+01  2.9744E+01  1.9148E+01  2.8554E+01  2.2063E+01
 EPSSHRINKSD(%)  2.6489E+01
 EPSSHRINKVR(%)  4.5962E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          700
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1286.5139464865417     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2853.4113023942632     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1566.8973559077215     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.75
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2853.411       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.07E+00  1.08E+00  1.12E+00  1.07E+00  1.11E+00  1.40E+00  1.87E+00  1.30E+00  1.08E+00  1.13E+00  1.04E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       89.853
Stop Time:
Sat Oct 23 20:35:51 CDT 2021
