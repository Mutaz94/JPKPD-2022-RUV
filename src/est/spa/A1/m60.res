Wed Sep 29 12:14:40 CDT 2021
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
$DATA ../../../../data/spa/A1/dat60.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m60.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -948.852162235675        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9209E+02 -2.4009E+01  2.4866E+01 -5.8374E+00  8.3356E+01  2.8672E+01 -2.1983E+01 -3.7261E+01 -1.0345E+01 -1.0547E+01
            -1.3104E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1392.19729408546        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0462E+00  1.0977E+00  9.5894E-01  1.0081E+00  9.6979E-01  1.1295E+00  1.0716E+00  1.1147E+00  9.4644E-01  7.9607E-01
             2.1407E+00
 PARAMETER:  1.4515E-01  1.9321E-01  5.8071E-02  1.0804E-01  6.9320E-02  2.2180E-01  1.6915E-01  2.0855E-01  4.4950E-02 -1.2806E-01
             8.6113E-01
 GRADIENT:   1.8094E+02  1.2299E+01  1.0087E+01  1.1304E+01  5.0832E+00  5.5240E+01 -5.8754E+00 -1.0083E+01 -8.0488E-01  7.6641E+00
            -1.0981E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1400.20236401556        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0467E+00  9.2579E-01  5.6485E-01  1.0755E+00  6.5478E-01  9.9492E-01  1.4751E+00  8.1250E-01  8.9451E-01  2.6678E-01
             2.1608E+00
 PARAMETER:  1.4566E-01  2.2894E-02 -4.7119E-01  1.7281E-01 -3.2346E-01  9.4907E-02  4.8872E-01 -1.0764E-01 -1.1485E-02 -1.2213E+00
             8.7046E-01
 GRADIENT:   1.8169E+02  2.8558E+01  1.4975E+01  4.7079E+01 -1.1035E+01  4.8219E-03  1.9537E+01 -4.3472E+00  1.9444E+01  1.2166E+00
            -9.1868E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1411.80106304925        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      333
 NPARAMETR:  1.0163E+00  8.9546E-01  8.5873E-01  1.1178E+00  8.2110E-01  9.9142E-01  1.3937E+00  1.1857E+00  7.4930E-01  2.1053E-01
             2.5947E+00
 PARAMETER:  1.1617E-01 -1.0415E-02 -5.2297E-02  2.1134E-01 -9.7113E-02  9.1384E-02  4.3194E-01  2.7032E-01 -1.8861E-01 -1.4581E+00
             1.0535E+00
 GRADIENT:  -2.4250E+00  3.7632E-01  2.9339E+00 -1.9608E+00 -7.6533E+00  1.7511E-01 -6.7509E-01  8.0501E-01 -1.4298E+00  9.1658E-01
             4.7662E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1413.14250003227        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      510
 NPARAMETR:  1.0184E+00  1.2092E+00  6.1491E-01  9.0955E-01  8.3479E-01  1.0001E+00  1.0895E+00  1.1507E+00  8.3892E-01  5.4680E-02
             2.5757E+00
 PARAMETER:  1.1819E-01  2.8998E-01 -3.8627E-01  5.1963E-03 -8.0580E-02  1.0006E-01  1.8574E-01  2.4038E-01 -7.5641E-02 -2.8063E+00
             1.0461E+00
 GRADIENT:  -3.5473E+00 -9.1535E+00 -4.2373E+00  4.2057E+00  1.3041E+01  2.1711E+00 -2.1919E-01 -2.5026E-02 -8.0116E-01  7.0872E-02
             1.6353E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1417.76979773998        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      694
 NPARAMETR:  1.0181E+00  1.3506E+00  2.2918E-01  7.3147E-01  6.2844E-01  9.8314E-01  8.8631E-01  7.7310E-01  8.4549E-01  1.1512E-02
             2.4814E+00
 PARAMETER:  1.1791E-01  4.0056E-01 -1.3732E+00 -2.1270E-01 -3.6452E-01  8.2997E-02 -2.0693E-02 -1.5735E-01 -6.7836E-02 -4.3643E+00
             1.0088E+00
 GRADIENT:   1.5256E+01 -1.6001E+00 -4.4262E-01  2.5875E+01  6.6877E+00 -4.1620E+00 -7.0874E+00 -2.9924E+00 -4.4244E-01  6.3707E-03
            -2.6449E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1464.92536588685        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      874
 NPARAMETR:  1.0079E+00  1.6275E+00  1.8243E-01  5.2192E-01  7.5595E-01  9.7046E-01  8.1578E-01  2.4315E+00  9.5990E-01  1.0000E-02
             2.1514E+00
 PARAMETER:  1.0786E-01  5.8702E-01 -1.6014E+00 -5.5025E-01 -1.7978E-01  7.0019E-02 -1.0362E-01  9.8852E-01  5.9075E-02 -5.5945E+00
             8.6611E-01
 GRADIENT:   1.8773E+01 -5.5147E+01  1.4059E+01 -3.7483E+01 -6.8044E+00  1.2852E-01 -1.1535E+01 -3.2447E+00 -4.6230E+00  0.0000E+00
             4.7188E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1469.88681195597        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1050
 NPARAMETR:  9.9841E-01  1.6617E+00  1.6585E-01  5.1825E-01  7.5103E-01  9.7345E-01  8.3538E-01  2.3942E+00  1.1285E+00  1.0000E-02
             1.9108E+00
 PARAMETER:  9.8408E-02  6.0784E-01 -1.6967E+00 -5.5730E-01 -1.8630E-01  7.3092E-02 -7.9864E-02  9.7305E-01  2.2088E-01 -5.3203E+00
             7.4751E-01
 GRADIENT:   4.3615E+00 -4.8929E+00  5.5486E+00 -1.1412E+01 -9.2018E+00 -2.6708E-01 -3.0774E+00  4.9520E-01 -4.1765E-01  0.0000E+00
             8.7668E-01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1470.04418346852        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1177
 NPARAMETR:  9.9585E-01  1.6645E+00  1.5831E-01  5.2137E-01  7.4723E-01  9.7484E-01  8.4376E-01  2.3449E+00  1.1324E+00  1.0000E-02
             1.9039E+00
 PARAMETER:  9.5845E-02  6.0953E-01 -1.7432E+00 -5.5129E-01 -1.9138E-01  7.4519E-02 -6.9890E-02  9.5226E-01  2.2433E-01 -5.3451E+00
             7.4393E-01
 GRADIENT:  -6.6282E-01  1.7402E+00  1.1921E-01  1.1690E+00 -3.7593E-01  1.0235E-02 -1.5540E-01 -5.4622E-02 -5.5360E-02  0.0000E+00
            -6.7855E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1177
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.1575E-04 -1.2546E-02 -8.5029E-03  1.1599E-02 -5.0433E-04
 SE:             2.9480E-02  2.6693E-02  1.7697E-02  2.0587E-02  2.7233E-04
 N:                     100         100         100         100         100

 P VAL.:         9.7792E-01  6.3836E-01  6.3090E-01  5.7315E-01  6.4039E-02

 ETASHRINKSD(%)  1.2377E+00  1.0576E+01  4.0712E+01  3.1030E+01  9.9088E+01
 ETASHRINKVR(%)  2.4600E+00  2.0034E+01  6.4849E+01  5.2432E+01  9.9992E+01
 EBVSHRINKSD(%)  1.4429E+00  1.0690E+01  4.1064E+01  3.1415E+01  9.9086E+01
 EBVSHRINKVR(%)  2.8650E+00  2.0237E+01  6.5265E+01  5.2960E+01  9.9992E+01
 RELATIVEINF(%)  9.5682E+01  1.2893E+01  1.3361E+01  7.7877E+00  1.3228E-03
 EPSSHRINKSD(%)  3.9264E+01
 EPSSHRINKVR(%)  6.3111E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1470.0441834685221     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -734.89335690478390     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.45
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.51
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1470.044       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.96E-01  1.66E+00  1.58E-01  5.21E-01  7.47E-01  9.75E-01  8.44E-01  2.34E+00  1.13E+00  1.00E-02  1.90E+00
 


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
+        1.15E+03
 
 TH 2
+        1.58E+01  3.26E+04
 
 TH 3
+       -5.86E+01  2.87E+03  2.81E+03
 
 TH 4
+       -5.35E+01 -4.04E+02 -1.34E+03  1.49E+03
 
 TH 5
+       -5.18E+01 -5.33E+02 -1.32E+03  7.46E+02  1.63E+06
 
 TH 6
+        6.02E-01  2.19E+01  3.88E+00 -1.25E+01 -3.76E+00  1.96E+02
 
 TH 7
+        2.57E+00  6.02E+01 -3.16E+01 -1.78E+01  1.59E+01  1.81E+00  1.82E+02
 
 TH 8
+        1.79E+01  1.44E+04  1.11E+03 -2.80E+02 -3.13E+01  1.37E+01  2.17E+01  6.46E+03
 
 TH 9
+        1.34E+00 -9.55E+00  8.99E+00  2.97E+01 -9.16E+00 -5.79E-01  1.50E+01  2.91E+00  3.45E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.36E+05 -1.26E+02  8.43E+04  4.58E+02  1.64E+05 -1.77E+01 -2.66E+01 -1.04E+04  1.32E+01  0.00E+00  1.67E+04
 
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
 #CPUT: Total CPU Time in Seconds,       21.027
Stop Time:
Wed Sep 29 12:15:02 CDT 2021
