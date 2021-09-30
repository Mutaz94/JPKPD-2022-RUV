Wed Sep 29 11:33:05 CDT 2021
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
$DATA ../../../../data/spa/B/dat80.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m80.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1716.86235666759        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8700E+02  6.6501E+00  1.1245E+01  1.9482E+01 -3.1855E+01  8.0039E+01  1.5254E+01  6.8688E+00  4.0093E+01  2.0171E+01
            -3.1954E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1726.52110764010        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.0265E+00  1.0568E+00  1.0294E+00  1.0164E+00  1.0580E+00  8.8102E-01  9.3546E-01  9.6259E-01  8.4090E-01  9.1070E-01
             1.1286E+00
 PARAMETER:  1.2615E-01  1.5527E-01  1.2895E-01  1.1630E-01  1.5638E-01 -2.6680E-02  3.3281E-02  6.1871E-02 -7.3284E-02  6.4631E-03
             2.2098E-01
 GRADIENT:  -1.0279E+00  1.4872E+01  1.9788E+00  1.9277E+01  3.2057E+00 -1.2349E+01  2.1059E+00  6.0594E-01  1.3709E-01  2.4121E+00
             1.1417E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1727.36795942160        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  1.0317E+00  1.0464E+00  8.6935E-01  1.0059E+00  9.6091E-01  9.0943E-01  9.7385E-01  8.3631E-01  8.0955E-01  7.5729E-01
             1.1247E+00
 PARAMETER:  1.3125E-01  1.4536E-01 -4.0008E-02  1.0589E-01  6.0123E-02  5.0629E-03  7.3505E-02 -7.8752E-02 -1.1128E-01 -1.7800E-01
             2.1750E-01
 GRADIENT:   9.5733E+00  4.0281E+00 -1.6007E+00  9.0836E+00  6.1558E+00 -2.3415E-01  1.5444E-02  1.3063E+00 -3.2347E+00 -3.1957E+00
             1.0595E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1728.62219068566        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      547
 NPARAMETR:  1.0301E+00  1.2432E+00  6.1327E-01  8.6173E-01  9.0941E-01  9.1557E-01  8.7981E-01  3.3126E-01  8.9261E-01  7.2766E-01
             1.0846E+00
 PARAMETER:  1.2969E-01  3.1770E-01 -3.8895E-01 -4.8819E-02  5.0462E-03  1.1795E-02 -2.8046E-02 -1.0049E+00 -1.3604E-02 -2.1792E-01
             1.8120E-01
 GRADIENT:  -2.6745E-01  1.3093E+00 -1.5972E+00  5.0812E+00  3.5236E+00  9.2936E-01  6.9992E-01  3.1945E-01 -9.5304E-02 -7.5216E-01
            -1.8268E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1728.80983578663        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      724
 NPARAMETR:  1.0304E+00  1.3899E+00  5.2752E-01  7.6347E-01  9.3594E-01  9.1309E-01  8.0032E-01  1.3983E-01  9.6295E-01  7.2214E-01
             1.0905E+00
 PARAMETER:  1.2993E-01  4.2923E-01 -5.3956E-01 -1.6988E-01  3.3792E-02  9.0794E-03 -1.2275E-01 -1.8673E+00  6.2250E-02 -2.2553E-01
             1.8663E-01
 GRADIENT:  -6.1248E-01  8.7909E-01  1.2722E-01  7.0482E-01 -5.9851E-01 -4.1044E-01 -6.9338E-01  6.7894E-02 -5.3521E-01 -2.0491E-01
             8.0073E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1728.83730486493        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  1.0315E+00  1.4017E+00  5.2226E-01  7.5439E-01  9.4102E-01  9.1456E-01  7.9987E-01  4.8373E-02  9.7563E-01  7.2565E-01
             1.0904E+00
 PARAMETER:  1.3101E-01  4.3766E-01 -5.4959E-01 -1.8184E-01  3.9214E-02  1.0687E-02 -1.2331E-01 -2.9288E+00  7.5328E-02 -2.2069E-01
             1.8655E-01
 GRADIENT:   2.3493E+00 -2.2540E+00 -3.6284E-01 -1.1679E+00  4.4457E-01  2.1486E-01  4.6568E-01  8.5940E-03  4.6204E-01  1.3810E-01
             2.5223E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1728.83985030977        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1061
 NPARAMETR:  1.0305E+00  1.4026E+00  5.2238E-01  7.5488E-01  9.4112E-01  9.1399E-01  7.9753E-01  2.9139E-02  9.7343E-01  7.2630E-01
             1.0901E+00
 PARAMETER:  1.3004E-01  4.3836E-01 -5.4936E-01 -1.8119E-01  3.9318E-02  1.0069E-02 -1.2623E-01 -3.4357E+00  7.3071E-02 -2.1980E-01
             1.8623E-01
 GRADIENT:  -2.9415E-01 -8.4240E-02 -8.6620E-02  1.9630E-01  2.6942E-02 -3.5025E-02 -4.3309E-02  3.0230E-03 -3.1180E-02 -1.4895E-02
            -1.9882E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1728.84082980988        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1237
 NPARAMETR:  1.0306E+00  1.4037E+00  5.2295E-01  7.5433E-01  9.4228E-01  9.1417E-01  7.9723E-01  1.0000E-02  9.7453E-01  7.2759E-01
             1.0901E+00
 PARAMETER:  1.3012E-01  4.3910E-01 -5.4827E-01 -1.8192E-01  4.0546E-02  1.0264E-02 -1.2661E-01 -4.9876E+00  7.4201E-02 -2.1801E-01
             1.8630E-01
 GRADIENT:  -7.1370E-02 -6.2233E-02 -8.0091E-03  1.1624E-01 -8.3966E-03  4.2947E-02  6.0150E-03  0.0000E+00  2.6654E-03 -8.6269E-03
            -3.9911E-03

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1728.84175955129        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:     1367
 NPARAMETR:  1.0320E+00  1.4031E+00  5.2299E-01  7.5403E-01  9.4247E-01  9.1435E-01  7.9714E-01  1.0000E-02  9.7474E-01  7.2791E-01
             1.0902E+00
 PARAMETER:  1.3149E-01  4.3868E-01 -5.4819E-01 -1.8232E-01  4.0751E-02  1.0463E-02 -1.2673E-01 -5.1799E+00  7.4420E-02 -2.1758E-01
             1.8632E-01
 GRADIENT:   3.6051E+00 -1.5946E+00 -1.8159E-01 -7.6323E-01  4.8497E-01  1.2428E-01 -6.5523E-03  0.0000E+00  3.6516E-02  3.1124E-02
             2.2771E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1367
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.1824E-04 -1.5147E-02 -3.2714E-04  1.0393E-02 -2.3167E-02
 SE:             2.9808E-02  2.3437E-02  1.4868E-04  2.3730E-02  2.0696E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8881E-01  5.1811E-01  2.7782E-02  6.6140E-01  2.6297E-01

 ETASHRINKSD(%)  1.3849E-01  2.1483E+01  9.9502E+01  2.0501E+01  3.0666E+01
 ETASHRINKVR(%)  2.7680E-01  3.8350E+01  9.9998E+01  3.6800E+01  5.1927E+01
 EBVSHRINKSD(%)  5.7405E-01  2.1423E+01  9.9535E+01  2.0976E+01  3.0647E+01
 EBVSHRINKVR(%)  1.1448E+00  3.8256E+01  9.9998E+01  3.7552E+01  5.1901E+01
 RELATIVEINF(%)  9.8801E+01  2.7175E+00  1.3631E-04  2.9217E+00  4.4376E+00
 EPSSHRINKSD(%)  4.2787E+01
 EPSSHRINKVR(%)  6.7267E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1728.8417595512924     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -993.69093298755422     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.41
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.71
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1728.842       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.40E+00  5.23E-01  7.54E-01  9.42E-01  9.14E-01  7.97E-01  1.00E-02  9.75E-01  7.28E-01  1.09E+00
 


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
+        1.23E+03
 
 TH 2
+       -1.05E+01  5.64E+02
 
 TH 3
+        1.18E+01  3.34E+02  9.32E+02
 
 TH 4
+       -2.51E+01  3.94E+02 -5.71E+02  1.30E+03
 
 TH 5
+       -6.27E+00 -4.35E+02 -8.92E+02  4.77E+02  1.13E+03
 
 TH 6
+       -1.81E+00 -1.60E+00  2.94E+00 -4.81E+00 -1.09E+00  2.33E+02
 
 TH 7
+        5.91E-01  2.29E+01 -3.02E+01 -1.24E+01 -2.07E+00  3.64E-01  1.26E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.67E+00 -2.38E+01 -4.57E+01  5.92E+01 -7.64E+00 -4.78E-01  1.99E+01  0.00E+00  9.15E+01
 
 TH10
+       -4.35E-01 -1.30E+01 -5.52E+01 -2.15E+01 -7.00E+01 -4.98E-02  3.07E+01  0.00E+00  1.34E+01  9.30E+01
 
 TH11
+       -7.76E+00 -1.97E+01 -3.28E+01 -3.25E+00 -8.78E+00  2.69E+00  1.02E+01  0.00E+00  1.29E+01  2.50E+01  1.82E+02
 
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
 #CPUT: Total CPU Time in Seconds,       23.188
Stop Time:
Wed Sep 29 11:33:30 CDT 2021
