Sat Sep 18 10:48:59 CDT 2021
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
$DATA ../../../../data/spa/A3/dat95.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       18 SEP 2021
Days until program expires : 211
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
 NO. OF SIG. FIGURES REQUIRED:            3
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
 RAW OUTPUT FILE (FILE): m95.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   68.4054858899652        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.5734E+01  1.7550E+01  1.7518E+01 -7.7394E+01  2.4926E+02  2.8877E+01 -1.1976E+02 -4.5633E+00 -2.7409E+02 -1.5995E+02
            -2.9079E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1217.07853712010        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0176E+00  8.7557E-01  1.1147E+00  1.1059E+00  8.8490E-01  7.3183E-01  1.0621E+00  8.8073E-01  1.5641E+00  7.5349E-01
             4.4699E+00
 PARAMETER:  1.1742E-01 -3.2880E-02  2.0855E-01  2.0069E-01 -2.2280E-02 -2.1221E-01  1.6025E-01 -2.7003E-02  5.4734E-01 -1.8303E-01
             1.5974E+00
 GRADIENT:  -8.0202E+00 -2.7676E+01 -1.7154E+00 -4.1233E+01  4.8898E+00 -3.5076E+01  9.9473E+00  6.6391E+00  2.9216E+01  1.6724E+01
             9.2083E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1243.21462215972        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0221E+00  4.2774E-01  3.7982E-01  1.4150E+00  3.6793E-01  8.6353E-01  1.1632E+00  1.9875E-01  1.2820E+00  1.9848E-01
             4.4195E+00
 PARAMETER:  1.2190E-01 -7.4923E-01 -8.6805E-01  4.4715E-01 -8.9985E-01 -4.6731E-02  2.5116E-01 -1.5157E+00  3.4840E-01 -1.5170E+00
             1.5860E+00
 GRADIENT:  -4.4563E+01  2.1712E+01 -4.3285E+01  1.0076E+02  4.1638E+01 -3.7909E+00  3.3929E+00  1.2585E+00 -1.1249E+01  2.6248E+00
             1.3116E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1267.92021952859        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0111E+00  3.4945E-01  2.9875E-01  1.2559E+00  2.9622E-01  8.9998E-01  1.0069E+00  1.0018E-02  1.3360E+00  2.4458E-01
             3.4137E+00
 PARAMETER:  1.1103E-01 -9.5140E-01 -1.1081E+00  3.2784E-01 -1.1167E+00 -5.3872E-03  1.0686E-01 -4.5033E+00  3.8965E-01 -1.3082E+00
             1.3278E+00
 GRADIENT:  -1.8179E+01  2.2001E+01  1.4435E+01  2.6355E+01 -3.5668E+01  1.0542E+00 -1.2258E+00 -8.3015E-05 -1.6058E+01 -8.2860E-01
            -2.3308E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1274.23218358158        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.0113E+00  1.4963E-01  3.2377E-01  1.3198E+00  2.8287E-01  8.9017E-01  1.2768E+00  1.0000E-02  1.3190E+00  2.2514E-01
             3.5279E+00
 PARAMETER:  1.1126E-01 -1.7996E+00 -1.0277E+00  3.7749E-01 -1.1628E+00 -1.6341E-02  3.4434E-01 -5.4370E+00  3.7687E-01 -1.3911E+00
             1.3607E+00
 GRADIENT:   1.3307E+00  6.9537E+00  4.2514E+01  4.5138E+00 -6.8676E+01  2.1580E+00 -3.3195E-01  0.0000E+00 -1.7777E+00 -7.2224E-01
             6.8271E-03

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1277.30007975661        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0003E+00  1.3248E-02  3.2991E-01  1.3731E+00  2.8074E-01  8.8035E-01  1.7959E+00  1.0000E-02  1.2612E+00  2.0963E-01
             3.5176E+00
 PARAMETER:  1.0026E-01 -4.2239E+00 -1.0089E+00  4.1704E-01 -1.1703E+00 -2.7439E-02  6.8551E-01 -7.7394E+00  3.3205E-01 -1.4624E+00
             1.3578E+00
 GRADIENT:   2.0526E+00  1.5716E-01  3.1683E+00  1.4608E+00 -6.0697E+00  4.4773E-01  8.4077E-04  0.0000E+00 -3.3206E-01  1.2715E-01
             8.2602E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1277.38437183962        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      443
 NPARAMETR:  9.9900E-01  1.0000E-02  3.3215E-01  1.3753E+00  2.8241E-01  8.7859E-01  1.9591E+00  1.0000E-02  1.2567E+00  1.8653E-01
             3.5220E+00
 PARAMETER:  9.9001E-02 -4.8200E+00 -1.0022E+00  4.1868E-01 -1.1644E+00 -2.9441E-02  7.7251E-01 -8.2175E+00  3.2845E-01 -1.5792E+00
             1.3590E+00
 GRADIENT:   4.1136E-01  0.0000E+00  5.9628E-01  3.2022E-01 -1.0717E+00  5.7254E-02  5.3525E-04  0.0000E+00 -4.1534E-01  5.1481E-02
            -1.1104E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1277.49879592275        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      580
 NPARAMETR:  9.9959E-01  1.0000E-02  3.4521E-01  1.3935E+00  2.9045E-01  8.7681E-01  2.4207E+00  1.0000E-02  1.2489E+00  1.3540E-01
             3.5458E+00
 PARAMETER:  9.9592E-02 -5.9694E+00 -9.6360E-01  4.3179E-01 -1.1363E+00 -3.1463E-02  9.8406E-01 -8.9684E+00  3.2228E-01 -1.8995E+00
             1.3658E+00
 GRADIENT:  -1.2314E-03  0.0000E+00  8.0709E-03  3.9766E-02 -1.5021E-02  2.1231E-03  5.7185E-04  0.0000E+00  7.8580E-03 -5.8505E-04
             5.6425E-02

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1277.49880486622        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:      706
 NPARAMETR:  9.9959E-01  1.0000E-02  3.4517E-01  1.3934E+00  2.9043E-01  8.7684E-01  2.4133E+00  1.0000E-02  1.2489E+00  1.3606E-01
             3.5454E+00
 PARAMETER:  9.9589E-02 -5.9576E+00 -9.6372E-01  4.3173E-01 -1.1364E+00 -3.1434E-02  9.8099E-01 -8.9608E+00  3.2230E-01 -1.8946E+00
             1.3657E+00
 GRADIENT:  -1.3281E-03  0.0000E+00 -2.5256E-02 -1.0468E-02  3.8053E-02  5.9549E-03  5.7444E-04  0.0000E+00  3.0127E-03  1.0481E-05
             1.3847E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      706
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -8.1401E-04 -8.2288E-04  1.4903E-04 -1.3351E-02  1.4444E-03
 SE:             2.8229E-02  3.5567E-04  2.4022E-04  2.7082E-02  4.8409E-03
 N:                     100         100         100         100         100

 P VAL.:         9.7700E-01  2.0689E-02  5.3502E-01  6.2202E-01  7.6542E-01

 ETASHRINKSD(%)  5.4302E+00  9.8808E+01  9.9195E+01  9.2728E+00  8.3782E+01
 ETASHRINKVR(%)  1.0565E+01  9.9986E+01  9.9994E+01  1.7686E+01  9.7370E+01
 EBVSHRINKSD(%)  5.1733E+00  9.8872E+01  9.9153E+01  8.2071E+00  8.3748E+01
 EBVSHRINKVR(%)  1.0079E+01  9.9987E+01  9.9993E+01  1.5741E+01  9.7359E+01
 RELATIVEINF(%)  7.5498E+01  1.5538E-03  1.9307E-04  3.3589E+01  4.7748E-02
 EPSSHRINKSD(%)  2.2976E+01
 EPSSHRINKVR(%)  4.0673E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1277.4988048662190     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -542.34797830248078     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.04
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.25
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1277.499       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.00E-02  3.45E-01  1.39E+00  2.90E-01  8.77E-01  2.41E+00  1.00E-02  1.25E+00  1.36E-01  3.55E+00
 


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
+        1.30E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.54E+02  0.00E+00  6.11E+03
 
 TH 4
+       -4.27E+01  0.00E+00 -1.10E+02  3.08E+02
 
 TH 5
+        4.01E+02  0.00E+00 -1.01E+04 -3.65E+02  1.94E+04
 
 TH 6
+       -1.65E+01  0.00E+00  3.84E+01 -1.25E+01 -1.24E+01  1.93E+02
 
 TH 7
+       -1.52E-01  0.00E+00  9.58E-02  1.28E-02  2.97E-01 -2.11E-01  1.43E-02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.24E+01  0.00E+00  1.88E+01 -8.38E+00  1.14E+02  2.89E+00 -4.19E-02  0.00E+00  8.95E+01
 
 TH10
+       -1.05E+01  0.00E+00 -6.78E+01 -6.02E+00  1.47E+02  3.55E+00  2.99E-01  0.00E+00  9.42E-01  9.35E+00
 
 TH11
+       -2.25E+01  0.00E+00 -2.05E+01 -4.34E+00  4.52E+01  4.72E+00  9.62E-03  0.00E+00  4.15E+00  1.02E+01  3.35E+01
 
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
 #CPUT: Total CPU Time in Seconds,       13.364
Stop Time:
Sat Sep 18 10:49:14 CDT 2021
