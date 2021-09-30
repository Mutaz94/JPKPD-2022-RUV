Wed Sep 29 23:01:31 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat7.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m7.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -800.662166897200        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8598E+02  6.3856E+01  1.2422E+02  1.0560E+02  1.7085E+02  5.3257E+00 -2.1013E+01 -2.6977E+02 -2.9859E+01 -8.1991E+01
            -2.1096E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1630.23626377497        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:       90
 NPARAMETR:  9.7608E-01  1.0331E+00  9.9250E-01  9.8510E-01  1.0006E+00  1.0799E+00  8.5842E-01  1.0252E+00  7.4619E-01  8.2109E-01
             2.9759E+00
 PARAMETER:  7.5791E-02  1.3260E-01  9.2468E-02  8.4985E-02  1.0056E-01  1.7683E-01 -5.2657E-02  1.2484E-01 -1.9278E-01 -9.7119E-02
             1.1905E+00
 GRADIENT:   1.3645E+01 -2.6978E+01 -1.5345E+01 -2.0784E+01  2.1787E+01  1.5140E+01 -3.5059E-01  2.8456E+00 -3.2424E+00  1.3488E+01
             5.3926E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1634.49794417253        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  9.9641E-01  9.7018E-01  6.1627E-01  9.9907E-01  7.0365E-01  1.0769E+00  7.9334E-01  8.7305E-01  8.3476E-01  4.0110E-01
             2.9476E+00
 PARAMETER:  9.6399E-02  6.9725E-02 -3.8407E-01  9.9070E-02 -2.5147E-01  1.7409E-01 -1.3150E-01 -3.5763E-02 -8.0614E-02 -8.1354E-01
             1.1810E+00
 GRADIENT:   4.1591E+01  2.0857E+01  2.4470E+01  5.7834E+00 -4.2704E+01  1.0382E+01 -2.8271E+00  1.8188E+00  3.5393E+00  3.4477E+00
             5.6505E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1637.69578541512        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      320
 NPARAMETR:  9.9684E-01  9.4517E-01  4.4005E-01  9.8042E-01  5.9925E-01  1.0824E+00  9.5990E-01  6.9448E-01  8.0194E-01  3.3299E-01
             2.8115E+00
 PARAMETER:  9.6831E-02  4.3611E-02 -7.2088E-01  8.0224E-02 -4.1207E-01  1.7921E-01  5.9070E-02 -2.6459E-01 -1.2072E-01 -9.9965E-01
             1.1337E+00
 GRADIENT:  -5.5868E+00  1.5029E+00 -9.9605E-01  8.1167E+00 -2.6377E+00  6.0174E-01  1.0760E+00  2.9667E-01  2.8245E+00  2.7027E+00
             1.6821E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1639.69863288303        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      496
 NPARAMETR:  9.9784E-01  1.1844E+00  4.3981E-01  8.6169E-01  7.1335E-01  1.0816E+00  8.2337E-01  9.2362E-01  8.4252E-01  1.8827E-01
             2.7822E+00
 PARAMETER:  9.7833E-02  2.6923E-01 -7.2142E-01 -4.8855E-02 -2.3779E-01  1.7848E-01 -9.4347E-02  2.0546E-02 -7.1353E-02 -1.5699E+00
             1.1232E+00
 GRADIENT:  -4.1211E+00  7.9769E+00 -2.2015E+00  2.0242E+01  3.7799E+00 -4.5871E-01 -6.9689E-01 -1.2182E+00 -1.3108E+00  9.5388E-01
            -4.8198E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1646.51407834114        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      676
 NPARAMETR:  1.0017E+00  1.7082E+00  3.7995E-01  5.8141E-01  9.5328E-01  1.0896E+00  6.5101E-01  2.1710E+00  1.0521E+00  1.5897E-02
             2.7473E+00
 PARAMETER:  1.0171E-01  6.3546E-01 -8.6771E-01 -4.4229E-01  5.2152E-02  1.8581E-01 -3.2923E-01  8.7518E-01  1.5081E-01 -4.0417E+00
             1.1106E+00
 GRADIENT:  -1.5381E-01  8.6467E+01  6.4008E+00  4.1227E+01 -3.2260E+01  3.5493E-01 -3.3545E+00 -3.9814E+00 -5.1251E+00  6.8060E-03
             5.5035E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1652.44730875342        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      853
 NPARAMETR:  9.9623E-01  1.9395E+00  4.1374E-01  4.1434E-01  1.2001E+00  1.0798E+00  5.7789E-01  2.8141E+00  1.5068E+00  1.0000E-02
             2.7556E+00
 PARAMETER:  9.6227E-02  7.6242E-01 -7.8252E-01 -7.8107E-01  2.8238E-01  1.7678E-01 -4.4837E-01  1.1346E+00  5.1000E-01 -6.5783E+00
             1.1136E+00
 GRADIENT:  -4.8317E+00  6.2752E+00  3.6554E-02  8.9477E+00  6.1597E+00 -1.6368E+00 -7.3547E-01 -4.1988E+00 -1.1066E+00  0.0000E+00
            -1.2937E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1656.67583398129        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1044
 NPARAMETR:  9.9883E-01  2.0582E+00  3.9411E-01  3.3041E-01  1.2831E+00  1.0866E+00  5.5934E-01  3.3314E+00  1.8723E+00  1.0000E-02
             2.7561E+00
 PARAMETER:  9.8828E-02  8.2182E-01 -8.3112E-01 -1.0074E+00  3.4925E-01  1.8307E-01 -4.8099E-01  1.3034E+00  7.2717E-01 -7.9636E+00
             1.1138E+00
 GRADIENT:   7.4507E-01 -6.1453E+00  1.0299E+00  9.7379E-01  5.3388E+00  1.4236E-01 -3.1654E-01 -5.4935E+00  4.0487E-02  0.0000E+00
             9.9311E+00

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1656.70508749155        NO. OF FUNC. EVALS.:  64
 CUMULATIVE NO. OF FUNC. EVALS.:     1108
 NPARAMETR:  9.9861E-01  2.0715E+00  3.9543E-01  3.2613E-01  1.2829E+00  1.0861E+00  5.6282E-01  3.2999E+00  1.8767E+00  1.0000E-02
             2.7720E+00
 PARAMETER:  9.8391E-02  8.2249E-01 -8.3596E-01 -1.0103E+00  3.4576E-01  1.8266E-01 -4.7909E-01  1.3038E+00  7.2659E-01 -7.9636E+00
             1.1133E+00
 GRADIENT:  -1.9341E-01 -3.3218E+01 -3.2169E+01  2.4513E+01 -7.7177E+01  1.4823E-02 -2.2904E-01  1.8754E+01 -9.8785E-02  0.0000E+00
            -1.8897E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1108
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.2931E-03 -3.4147E-02 -2.6481E-02  4.0305E-02 -2.6448E-04
 SE:             2.9464E-02  2.0889E-02  1.4386E-02  1.9929E-02  1.5384E-04
 N:                     100         100         100         100         100

 P VAL.:         9.1101E-01  1.0211E-01  6.5650E-02  4.3136E-02  8.5580E-02

 ETASHRINKSD(%)  1.2920E+00  3.0019E+01  5.1806E+01  3.3235E+01  9.9485E+01
 ETASHRINKVR(%)  2.5673E+00  5.1026E+01  7.6773E+01  5.5424E+01  9.9997E+01
 EBVSHRINKSD(%)  1.8457E+00  2.7737E+01  5.7854E+01  3.9249E+01  9.9386E+01
 EBVSHRINKVR(%)  3.6574E+00  4.7780E+01  8.2237E+01  6.3094E+01  9.9996E+01
 RELATIVEINF(%)  9.6069E+01  7.8167E+00  6.1841E+00  6.2075E+00  1.2092E-03
 EPSSHRINKSD(%)  2.3557E+01
 EPSSHRINKVR(%)  4.1564E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1656.7050874915460     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -737.76655428687332     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.12
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.60
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1656.705       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.98E-01  2.06E+00  3.92E-01  3.29E-01  1.28E+00  1.09E+00  5.60E-01  3.33E+00  1.87E+00  1.00E-02  2.75E+00
 


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
+        9.13E+02
 
 TH 2
+       -7.28E+01  6.99E+02
 
 TH 3
+       -9.48E+01  1.72E+02  6.86E+03
 
 TH 4
+        1.86E+04 -7.26E+02 -2.44E+02  6.99E+03
 
 TH 5
+       -8.22E+01  7.92E+02  4.82E+03 -6.36E+01  3.86E+03
 
 TH 6
+        3.62E+00 -3.55E+01 -6.61E+01  9.42E+03 -6.27E+01  1.56E+02
 
 TH 7
+        6.83E+00 -5.79E+01 -4.96E+01  7.03E+03 -5.35E+01  3.68E+00  1.73E+02
 
 TH 8
+        4.81E-01 -1.64E+00 -1.02E+01 -6.13E+00  6.23E+00  4.82E-01 -1.52E+00  4.78E+01
 
 TH 9
+        8.92E-01 -1.50E+01 -2.62E+01  1.43E+03 -1.71E+01 -1.95E-02  2.08E+01  1.55E+00  9.87E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.22E+01 -2.10E+01 -2.85E+00  2.31E+01 -2.04E+01  1.94E+00  1.52E+01 -6.18E+01  4.52E+00  0.00E+00  1.55E+02
 
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
 #CPUT: Total CPU Time in Seconds,       27.796
Stop Time:
Wed Sep 29 23:02:00 CDT 2021
