Sat Sep 18 09:07:30 CDT 2021
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
$DATA ../../../../data/spa/A1/dat20.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m20.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1010.17122881357        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0483E+02 -5.3399E+00  2.2320E+01 -8.9332E+00  8.1734E+01 -2.2709E+00 -2.7823E+01 -1.9149E+01 -1.1455E+01 -5.6062E+01
            -1.1436E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1381.32488624321        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0256E+00  1.0317E+00  1.0105E+00  1.0512E+00  9.6668E-01  8.8083E-01  1.0066E+00  9.7502E-01  8.8112E-01  9.0307E-01
             3.1374E+00
 PARAMETER:  1.2529E-01  1.3125E-01  1.1044E-01  1.4993E-01  6.6114E-02 -2.6886E-02  1.0656E-01  7.4701E-02 -2.6562E-02 -1.9587E-03
             1.2434E+00
 GRADIENT:   8.6909E+01  1.7094E+01 -1.0786E+01  3.2527E+01 -5.9657E+00 -3.1202E+01  4.2884E+00  5.8002E+00  4.8979E+00  1.5948E+01
             9.1578E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1391.83372114512        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0171E+00  8.0709E-01  6.1887E-01  1.1657E+00  6.6933E-01  9.6157E-01  1.4065E+00  5.0233E-01  8.0488E-01  4.7807E-01
             2.8491E+00
 PARAMETER:  1.1700E-01 -1.1433E-01 -3.7986E-01  2.5328E-01 -3.0148E-01  6.0815E-02  4.4111E-01 -5.8850E-01 -1.1706E-01 -6.3800E-01
             1.1470E+00
 GRADIENT:   5.4641E+01  1.9854E+01 -4.3309E+01  8.5750E+01  4.7106E+01 -5.8224E+00  1.1092E+01  2.4485E+00  3.5731E+00  3.1277E+00
             5.5120E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1400.47434515027        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.8421E-01  6.7637E-01  4.8312E-01  1.1567E+00  5.2129E-01  9.8821E-01  1.4327E+00  1.7788E-01  7.8793E-01  5.4075E-01
             2.4222E+00
 PARAMETER:  8.4084E-02 -2.9102E-01 -6.2749E-01  2.4561E-01 -5.5144E-01  8.8139E-02  4.5957E-01 -1.6267E+00 -1.3834E-01 -5.1480E-01
             9.8467E-01
 GRADIENT:  -6.9004E+00  1.5587E+01 -6.0360E+00  3.1077E+01  1.8097E+00  4.5625E-01  2.2900E+00  9.1912E-02 -2.2974E+00 -3.9177E-01
            -4.8577E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1403.17326638200        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.8224E-01  4.0692E-01  5.4099E-01  1.2949E+00  4.8528E-01  9.7510E-01  2.0598E+00  4.0782E-02  7.4645E-01  6.0367E-01
             2.4421E+00
 PARAMETER:  8.2083E-02 -7.9915E-01 -5.1435E-01  3.5842E-01 -6.2303E-01  7.4790E-02  8.2261E-01 -3.0995E+00 -1.9243E-01 -4.0473E-01
             9.9286E-01
 GRADIENT:  -2.3967E+00  4.2515E+00 -8.7079E-01  1.0561E+01 -1.3156E+00 -8.7984E-01 -1.0621E+00  1.3935E-02 -2.9775E+00 -2.0708E-01
             4.3342E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1404.00205951928        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  9.7920E-01  2.5573E-01  5.4174E-01  1.3607E+00  4.5630E-01  9.7266E-01  2.7991E+00  1.0000E-02  7.4369E-01  6.3330E-01
             2.4272E+00
 PARAMETER:  7.8977E-02 -1.2636E+00 -5.1297E-01  4.0797E-01 -6.8461E-01  7.2284E-02  1.1293E+00 -4.5252E+00 -1.9614E-01 -3.5681E-01
             9.8676E-01
 GRADIENT:   4.7583E-01  1.1927E+00  8.3713E-01  3.1896E+00 -2.5724E+00 -1.6037E-01  1.2587E-01  0.0000E+00  1.8405E-01  1.4752E-01
            -5.7068E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1404.14301869094        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      443
 NPARAMETR:  9.7700E-01  2.0549E-01  5.4376E-01  1.3822E+00  4.4999E-01  9.7105E-01  3.1843E+00  1.0000E-02  7.3931E-01  6.3536E-01
             2.4304E+00
 PARAMETER:  7.6736E-02 -1.4824E+00 -5.0925E-01  4.2368E-01 -6.9854E-01  7.0622E-02  1.2582E+00 -5.1866E+00 -2.0204E-01 -3.5357E-01
             9.8804E-01
 GRADIENT:  -3.3331E-01  4.4787E-01  7.0478E-01  1.3456E+00 -1.3563E+00 -8.1921E-02  1.9842E-01  0.0000E+00 -1.2526E-01  1.5843E-01
             5.2356E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1404.35202238979        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      600
 NPARAMETR:  9.7856E-01  1.9530E-01  5.9316E-01  1.4077E+00  4.7872E-01  9.7009E-01  3.3658E+00  1.0000E-02  7.3503E-01  6.4554E-01
             2.4462E+00
 PARAMETER:  7.8322E-02 -1.5332E+00 -4.2229E-01  4.4194E-01 -6.3664E-01  6.9637E-02  1.3137E+00 -5.2600E+00 -2.0785E-01 -3.3766E-01
             9.9455E-01
 GRADIENT:  -6.3077E-02  2.1100E-03  8.6318E-02  4.3400E-02 -8.8213E-02 -1.2468E-02 -2.1993E-02  0.0000E+00  4.6025E-02  2.2554E-03
            -2.3835E-03

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1404.35204545077        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      657
 NPARAMETR:  9.7857E-01  1.9485E-01  5.9285E-01  1.4078E+00  4.7848E-01  9.7013E-01  3.3703E+00  1.0000E-02  7.3489E-01  6.4544E-01
             2.4462E+00
 PARAMETER:  7.8336E-02 -1.5355E+00 -4.2281E-01  4.4200E-01 -6.3713E-01  6.9670E-02  1.3150E+00 -5.2671E+00 -2.0803E-01 -3.3783E-01
             9.9454E-01
 GRADIENT:  -5.7329E-03  4.6568E-03  7.0411E-03 -2.0685E-02 -6.0685E-03 -2.3992E-03  8.8779E-03  0.0000E+00  1.9179E-03  4.9615E-04
            -1.2210E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      657
 NO. OF SIG. DIGITS IN FINAL EST.:  3.9
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7223E-03  3.9876E-02 -1.3579E-04 -3.0607E-02  2.9058E-03
 SE:             2.9162E-02  1.6549E-02  1.9450E-04  2.3856E-02  1.7800E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5290E-01  1.5971E-02  4.8509E-01  1.9949E-01  8.7032E-01

 ETASHRINKSD(%)  2.3022E+00  4.4559E+01  9.9348E+01  2.0080E+01  4.0368E+01
 ETASHRINKVR(%)  4.5513E+00  6.9263E+01  9.9996E+01  3.6128E+01  6.4440E+01
 EBVSHRINKSD(%)  2.3149E+00  5.5190E+01  9.9277E+01  1.6810E+01  3.7119E+01
 EBVSHRINKVR(%)  4.5763E+00  7.9921E+01  9.9995E+01  3.0795E+01  6.0459E+01
 RELATIVEINF(%)  9.3934E+01  6.3033E+00  1.9601E-04  2.1258E+01  1.4951E+00
 EPSSHRINKSD(%)  3.2597E+01
 EPSSHRINKVR(%)  5.4568E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1404.3520454507670     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -669.20121888702886     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.26
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.76
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1404.352       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.79E-01  1.95E-01  5.93E-01  1.41E+00  4.78E-01  9.70E-01  3.37E+00  1.00E-02  7.35E-01  6.45E-01  2.45E+00
 


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
+       -3.75E+01  1.01E+03
 
 TH 3
+        8.05E+00  9.74E+01  1.85E+03
 
 TH 4
+       -5.46E+01  1.68E+02 -1.39E+02  7.92E+02
 
 TH 5
+        6.33E+01 -4.12E+02 -2.90E+03 -1.52E+02  4.95E+03
 
 TH 6
+       -6.43E+00  3.12E+00  9.36E+00 -1.82E+01  4.46E+00  1.90E+02
 
 TH 7
+        2.95E+00  8.05E+01 -2.60E+01 -1.95E+01  3.81E+01  1.31E+00  9.20E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.54E+00 -5.29E+01  3.93E+01 -1.15E+01 -1.40E+01  6.50E+00 -3.58E+00  0.00E+00  2.16E+02
 
 TH10
+       -1.17E+01 -7.05E+01 -2.43E+00  8.81E+00 -4.21E+01 -4.48E+00 -7.17E+00  0.00E+00 -6.77E+00  1.02E+02
 
 TH11
+       -1.40E+01 -2.69E+01 -1.06E+01 -4.02E+00 -6.54E+00  3.12E+00 -2.11E+00  0.00E+00  1.02E+01  3.00E+01  5.07E+01
 
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
 #CPUT: Total CPU Time in Seconds,       13.084
Stop Time:
Sat Sep 18 09:07:45 CDT 2021
