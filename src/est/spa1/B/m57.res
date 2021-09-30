Wed Sep 29 21:14:36 CDT 2021
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
$DATA ../../../../data/spa1/B/dat57.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m57.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2150.49294511943        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6630E+02 -5.1594E+01 -4.9284E+01  1.2081E+01  6.9989E+01  6.9277E+01  1.1394E+01  1.2800E+01  3.5277E+01  2.3701E+01
            -7.2503E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2166.29533008696        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      196
 NPARAMETR:  1.0423E+00  1.0831E+00  1.0874E+00  1.0217E+00  1.0232E+00  8.9855E-01  9.3979E-01  9.4903E-01  8.7717E-01  8.7590E-01
             1.1045E+00
 PARAMETER:  1.4148E-01  1.7979E-01  1.8378E-01  1.2150E-01  1.2290E-01 -6.9758E-03  3.7903E-02  4.7688E-02 -3.1056E-02 -3.2505E-02
             1.9943E-01
 GRADIENT:   2.2274E+01  4.0475E+00 -9.4994E+00  2.0092E+01  2.8556E+01 -8.7921E+00  2.5246E+00  1.6678E+00 -6.1162E-02  2.3616E+00
             5.5193E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2167.61977106089        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  1.0461E+00  1.0695E+00  8.7521E-01  1.0132E+00  9.0831E-01  9.0894E-01  1.0222E+00  6.9726E-01  8.6924E-01  7.2995E-01
             1.1008E+00
 PARAMETER:  1.4505E-01  1.6717E-01 -3.3292E-02  1.1310E-01  3.8343E-03  4.5272E-03  1.2198E-01 -2.6059E-01 -4.0130E-02 -2.1478E-01
             1.9602E-01
 GRADIENT:   2.8776E+01  2.5148E+00 -9.1851E+00  1.8180E+01  1.9663E+01 -4.6754E+00  3.8546E+00  1.8827E+00  4.1774E+00 -2.7691E+00
             2.6083E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2168.63040724710        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      549
 NPARAMETR:  1.0349E+00  1.1899E+00  7.3009E-01  9.2341E-01  8.8320E-01  9.2255E-01  9.0471E-01  3.7445E-01  9.1086E-01  7.3852E-01
             1.0962E+00
 PARAMETER:  1.3430E-01  2.7384E-01 -2.1459E-01  2.0321E-02 -2.4206E-02  1.9383E-02 -1.3992E-04 -8.8231E-01  6.6319E-03 -2.0310E-01
             1.9188E-01
 GRADIENT:  -4.0362E+00  3.0914E+00  7.6473E-01  1.5515E+00 -4.6699E+00  6.7095E-01 -9.8313E-01  4.6263E-01 -9.6825E-01  4.2263E-01
             4.0843E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2168.77128252127        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      726
 NPARAMETR:  1.0370E+00  1.2752E+00  6.8512E-01  8.6623E-01  9.0576E-01  9.2175E-01  8.6401E-01  1.6968E-01  9.6101E-01  7.4984E-01
             1.0959E+00
 PARAMETER:  1.3630E-01  3.4310E-01 -2.7817E-01 -4.3604E-02  1.0182E-03  1.8522E-02 -4.6170E-02 -1.6738E+00  6.0232E-02 -1.8790E-01
             1.9154E-01
 GRADIENT:  -1.1550E-01 -2.0015E+00 -9.1896E-01 -1.8458E+00  1.2020E+00  3.7981E-02  2.5881E-01  8.4803E-02  1.5382E-01  2.6728E-01
            -2.2454E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2168.82460505889        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      902
 NPARAMETR:  1.0369E+00  1.2355E+00  6.8803E-01  8.9104E-01  8.8618E-01  9.2139E-01  8.8949E-01  3.4011E-02  9.3856E-01  7.4086E-01
             1.0959E+00
 PARAMETER:  1.3622E-01  3.1149E-01 -2.7392E-01 -1.5370E-02 -2.0832E-02  1.8124E-02 -1.7110E-02 -3.2811E+00  3.6594E-02 -1.9994E-01
             1.9154E-01
 GRADIENT:   2.3515E-02  1.4907E-01 -2.4365E-01  2.5294E-01 -3.6389E-01  8.2583E-04  6.1319E-02  2.8760E-03  1.1405E-01  2.1112E-01
             6.4950E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2168.82691888306        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     1076
 NPARAMETR:  1.0377E+00  1.2358E+00  6.8832E-01  8.9098E-01  8.8646E-01  9.2164E-01  8.8956E-01  1.6591E-02  9.3829E-01  7.3996E-01
             1.0959E+00
 PARAMETER:  1.3701E-01  3.1170E-01 -2.7351E-01 -1.5437E-02 -2.0518E-02  1.8398E-02 -1.7030E-02 -3.9989E+00  3.6309E-02 -2.0116E-01
             1.9157E-01
 GRADIENT:   2.0906E+00  3.6662E-01  1.9609E-02  3.1795E-01 -2.9931E-01  1.1251E-01  1.4518E-02  6.5241E-04  1.4784E-03 -1.6640E-02
            -1.5028E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2168.82719545299        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1253
 NPARAMETR:  1.0379E+00  1.2357E+00  6.8839E-01  8.9054E-01  8.8683E-01  9.2155E-01  8.8919E-01  1.0000E-02  9.3851E-01  7.4030E-01
             1.0959E+00
 PARAMETER:  1.3718E-01  3.1161E-01 -2.7340E-01 -1.5932E-02 -2.0107E-02  1.8304E-02 -1.7448E-02 -5.1743E+00  3.6538E-02 -2.0070E-01
             1.9159E-01
 GRADIENT:   2.5480E+00 -6.1609E-01 -1.6708E-01 -3.9250E-01  2.2305E-01  7.8639E-02 -1.9320E-02  0.0000E+00  5.6726E-04 -4.0140E-03
             1.5753E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -2168.82719545299        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1275
 NPARAMETR:  1.0379E+00  1.2357E+00  6.8839E-01  8.9054E-01  8.8683E-01  9.2155E-01  8.8919E-01  1.0000E-02  9.3851E-01  7.4030E-01
             1.0959E+00
 PARAMETER:  1.3718E-01  3.1161E-01 -2.7340E-01 -1.5932E-02 -2.0107E-02  1.8304E-02 -1.7448E-02 -5.1743E+00  3.6538E-02 -2.0070E-01
             1.9159E-01
 GRADIENT:   2.5480E+00 -6.1609E-01 -1.6708E-01 -3.9250E-01  2.2305E-01  7.8639E-02 -1.9320E-02  0.0000E+00  5.6726E-04 -4.0140E-03
             1.5753E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1275
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.9377E-04 -1.0318E-02 -3.1965E-04  6.0850E-03 -1.6963E-02
 SE:             2.9834E-02  2.1910E-02  1.7209E-04  2.4764E-02  2.1920E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9214E-01  6.3769E-01  6.3245E-02  8.0590E-01  4.3901E-01

 ETASHRINKSD(%)  5.1195E-02  2.6598E+01  9.9423E+01  1.7038E+01  2.6566E+01
 ETASHRINKVR(%)  1.0236E-01  4.6122E+01  9.9997E+01  3.1172E+01  4.6075E+01
 EBVSHRINKSD(%)  4.5625E-01  2.6514E+01  9.9470E+01  1.7311E+01  2.6464E+01
 EBVSHRINKVR(%)  9.1043E-01  4.5999E+01  9.9997E+01  3.1625E+01  4.5924E+01
 RELATIVEINF(%)  9.8849E+01  3.2087E+00  3.3707E-04  4.9155E+00  5.4630E+00
 EPSSHRINKSD(%)  3.2526E+01
 EPSSHRINKVR(%)  5.4472E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2168.8271954529855     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1249.8886622483128     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.07
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.71
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2168.827       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.24E+00  6.88E-01  8.91E-01  8.87E-01  9.22E-01  8.89E-01  1.00E-02  9.39E-01  7.40E-01  1.10E+00
 


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
+        1.20E+03
 
 TH 2
+       -8.96E+00  5.14E+02
 
 TH 3
+        1.05E+01  2.79E+02  7.23E+02
 
 TH 4
+       -8.09E+00  3.91E+02 -3.70E+02  1.09E+03
 
 TH 5
+        1.12E+00 -4.61E+02 -8.42E+02  4.13E+02  1.39E+03
 
 TH 6
+        8.23E-01 -1.85E+00  3.61E+00 -2.78E+00 -5.01E-01  2.30E+02
 
 TH 7
+        1.33E+00  2.83E+01 -9.59E+00 -9.12E+00 -9.07E+00  8.10E-01  7.78E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.16E+00 -2.66E+01 -2.71E+01  4.43E+01 -7.19E-01 -5.81E-01  2.47E+01  0.00E+00  1.15E+02
 
 TH10
+        8.96E-01 -9.72E+00 -6.42E+01 -1.58E+01 -5.99E+01  5.36E-01  2.57E+01  0.00E+00  9.76E+00  1.12E+02
 
 TH11
+       -9.44E+00 -1.74E+01 -2.54E+01 -9.26E+00  7.94E-01  1.80E+00  7.16E+00  0.00E+00  7.66E+00  2.74E+01  3.45E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       25.847
Stop Time:
Wed Sep 29 21:15:03 CDT 2021
