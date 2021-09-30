Wed Sep 29 21:46:06 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat2.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1907.26160454413        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5267E+02 -6.5392E+00 -8.0848E+00  2.1591E+01 -7.2782E+00  3.8985E+01  6.5069E+00  1.6582E+01  4.3215E+01  1.5766E+01
            -4.5407E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1969.75326600515        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0114E+00  1.1345E+00  1.2002E+00  9.7632E-01  1.1660E+00  1.0216E+00  9.2078E-01  8.4604E-01  6.9377E-01  7.9609E-01
             1.8558E+00
 PARAMETER:  1.1130E-01  2.2620E-01  2.8248E-01  7.6038E-02  2.5359E-01  1.2133E-01  1.7468E-02 -6.7188E-02 -2.6562E-01 -1.2805E-01
             7.1832E-01
 GRADIENT:   1.8134E+02  2.3202E+01 -6.9505E+00  2.3598E+01  1.5239E+01  2.6590E+01 -1.0810E+01  2.3896E+00 -1.4371E+01 -5.8499E-02
             1.5799E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1975.16270444814        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.9146E-01  9.9311E-01  8.2675E-01  1.0478E+00  9.3484E-01  9.6967E-01  1.4091E+00  2.6826E-01  5.4901E-01  7.8547E-01
             1.7253E+00
 PARAMETER:  9.1422E-02  9.3081E-02 -9.0257E-02  1.4673E-01  3.2624E-02  6.9200E-02  4.4296E-01 -1.2158E+00 -4.9964E-01 -1.4147E-01
             6.4538E-01
 GRADIENT:   1.4262E+02  2.4335E+01 -5.1164E+01  9.5461E+01  6.4083E+01  6.3810E+00  2.6903E+01  8.0840E-01 -1.0795E+01  7.1656E+00
             1.2756E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1988.42686211680        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      249
 NPARAMETR:  9.2475E-01  1.0287E+00  8.4151E-01  9.9890E-01  9.3581E-01  9.5156E-01  1.0669E+00  1.5505E-01  7.6927E-01  8.4943E-01
             1.4809E+00
 PARAMETER:  2.1769E-02  1.2828E-01 -7.2557E-02  9.8902E-02  3.3657E-02  5.0343E-02  1.6474E-01 -1.7640E+00 -1.6232E-01 -6.3191E-02
             4.9264E-01
 GRADIENT:  -1.5231E+02 -1.7028E+01 -1.2954E+01 -2.4663E+01  6.9364E+00 -2.4199E+01 -1.2442E-02  3.1106E-01 -9.0266E-01  9.9163E+00
             2.4396E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1995.13123707092        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      426
 NPARAMETR:  9.8714E-01  8.6567E-01  9.0405E-01  1.1203E+00  8.8073E-01  9.7204E-01  1.2059E+00  1.0629E-01  7.2733E-01  7.8273E-01
             1.4324E+00
 PARAMETER:  8.7055E-02 -4.4256E-02 -8.6629E-04  2.1359E-01 -2.7008E-02  7.1639E-02  2.8720E-01 -2.1415E+00 -2.1837E-01 -1.4497E-01
             4.5933E-01
 GRADIENT:   9.6772E+00  1.1420E+01  2.4071E+00  1.5388E+01 -4.5598E+00 -5.7030E+00 -2.8244E+00  9.5384E-02 -2.8308E+00 -3.6742E-01
            -1.0458E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1996.09001129692        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      603
 NPARAMETR:  9.8203E-01  6.5896E-01  8.4420E-01  1.2224E+00  7.6622E-01  9.8454E-01  1.5494E+00  2.8163E-02  6.8947E-01  7.0034E-01
             1.4411E+00
 PARAMETER:  8.1867E-02 -3.1709E-01 -6.9368E-02  3.0082E-01 -1.6629E-01  8.4418E-02  5.3786E-01 -3.4698E+00 -2.7183E-01 -2.5619E-01
             4.6543E-01
 GRADIENT:   1.9729E+00  2.6065E+00  1.2391E+00  1.5490E+00 -1.8550E+00  1.9090E-01  4.8869E-01  1.0096E-02  6.1815E-04 -3.9852E-01
            -9.4589E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1996.13430504442        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      778
 NPARAMETR:  9.8005E-01  5.9402E-01  8.5348E-01  1.2573E+00  7.4981E-01  9.8294E-01  1.6609E+00  1.6415E-02  6.8078E-01  7.0658E-01
             1.4421E+00
 PARAMETER:  7.9845E-02 -4.2084E-01 -5.8429E-02  3.2900E-01 -1.8794E-01  8.2791E-02  6.0738E-01 -4.0095E+00 -2.8452E-01 -2.4732E-01
             4.6613E-01
 GRADIENT:  -2.7421E-01 -4.8745E-03  8.2017E-01 -1.8861E+00 -1.5321E+00 -3.2259E-02  7.3219E-02  3.4594E-03  1.6602E-01  7.7399E-02
             3.8441E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1996.13574174638        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      954
 NPARAMETR:  9.8038E-01  5.9206E-01  8.5763E-01  1.2594E+00  7.5231E-01  9.8308E-01  1.6636E+00  1.0000E-02  6.8009E-01  7.0949E-01
             1.4416E+00
 PARAMETER:  8.0184E-02 -4.2415E-01 -5.3587E-02  3.3060E-01 -1.8461E-01  8.2937E-02  6.0901E-01 -5.0007E+00 -2.8553E-01 -2.4322E-01
             4.6577E-01
 GRADIENT:   6.4725E-01 -2.8463E-01 -7.0408E-01 -4.0511E-01  7.8734E-01  3.4866E-02  2.0439E-02  0.0000E+00  4.8413E-02 -1.2834E-02
            -2.5269E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1996.13574174638        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:      982
 NPARAMETR:  9.8038E-01  5.9279E-01  8.5848E-01  1.2589E+00  7.5194E-01  9.8311E-01  1.6646E+00  1.0000E-02  6.8001E-01  7.0962E-01
             1.4417E+00
 PARAMETER:  8.0184E-02 -4.2415E-01 -5.3587E-02  3.3060E-01 -1.8461E-01  8.2937E-02  6.0901E-01 -5.0007E+00 -2.8553E-01 -2.4322E-01
             4.6577E-01
 GRADIENT:  -3.4343E-03 -3.2052E-01 -4.1941E-01  8.2190E-01  8.0778E-01 -8.2682E-03 -4.3986E-02  0.0000E+00  1.7096E-02 -1.4237E-02
            -1.6175E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      982
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.5072E-04  2.1203E-02 -2.9368E-04 -2.0392E-02 -4.2104E-03
 SE:             2.9738E-02  2.0061E-02  1.9553E-04  2.4073E-02  2.0762E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7450E-01  2.9054E-01  1.3311E-01  3.9695E-01  8.3930E-01

 ETASHRINKSD(%)  3.7343E-01  3.2793E+01  9.9345E+01  1.9351E+01  3.0443E+01
 ETASHRINKVR(%)  7.4546E-01  5.4832E+01  9.9996E+01  3.4958E+01  5.1619E+01
 EBVSHRINKSD(%)  7.0779E-01  3.5267E+01  9.9302E+01  1.8012E+01  2.8989E+01
 EBVSHRINKVR(%)  1.4106E+00  5.8096E+01  9.9995E+01  3.2780E+01  4.9575E+01
 RELATIVEINF(%)  9.7799E+01  4.0134E+00  3.2979E-04  7.5207E+00  3.5190E+00
 EPSSHRINKSD(%)  2.9726E+01
 EPSSHRINKVR(%)  5.0616E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1996.1357417463837     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1077.1972085417110     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.63
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.14
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1996.136       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.80E-01  5.92E-01  8.58E-01  1.26E+00  7.52E-01  9.83E-01  1.66E+00  1.00E-02  6.80E-01  7.09E-01  1.44E+00
 


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
+        1.18E+03
 
 TH 2
+       -2.07E+01  5.00E+02
 
 TH 3
+        9.10E+00  2.41E+02  7.71E+02
 
 TH 4
+       -7.23E+00  4.46E+02 -3.04E+02  1.03E+03
 
 TH 5
+        1.92E+00 -5.11E+02 -1.16E+03  2.94E+02  1.96E+03
 
 TH 6
+        1.19E+00 -4.26E+00  3.00E+00 -4.12E+00 -3.38E+00  2.00E+02
 
 TH 7
+        1.23E+00  3.62E+01  5.48E+00 -6.39E+00 -8.12E+00 -6.99E-02  1.94E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.35E+00 -1.48E+01 -1.94E+01 -7.52E+00  2.89E+01 -8.67E-01  1.43E+01  0.00E+00  2.13E+02
 
 TH10
+        5.78E-01  8.52E+00 -4.83E+01 -2.74E+01 -3.46E+01  4.95E-01  8.64E+00  0.00E+00  6.96E+00  1.02E+02
 
 TH11
+       -1.07E+01 -1.03E+01 -3.08E+01 -1.85E+01  9.43E+00  2.38E+00  4.18E+00  0.00E+00  1.24E+01  3.00E+01  2.09E+02
 
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
 #CPUT: Total CPU Time in Seconds,       20.841
Stop Time:
Wed Sep 29 21:46:28 CDT 2021
