Thu Sep 30 01:51:54 CDT 2021
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
$DATA ../../../../data/spa1/TD2/dat8.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m8.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2114.14429553083        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7067E+02  2.2246E+01  1.4021E+01  6.3642E+01  1.0589E+01  3.3517E+01 -6.0717E-01 -6.5742E+00  1.5949E+01  1.7003E+00
             3.1312E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2115.96031667655        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.8055E-01  1.0881E+00  8.7391E-01  9.6245E-01  9.5316E-01  1.0748E+00  1.0728E+00  1.0962E+00  9.1161E-01  9.7050E-01
             9.6838E-01
 PARAMETER:  8.0355E-02  1.8443E-01 -3.4783E-02  6.1729E-02  5.2031E-02  1.7213E-01  1.7026E-01  1.9188E-01  7.4575E-03  7.0060E-02
             6.7873E-02
 GRADIENT:  -7.7889E+00  1.5232E+01  3.7745E+00  1.2041E+01 -2.7093E+01  1.7347E+01 -1.8142E+00  2.7755E+00 -7.0678E+00  1.2299E+01
             7.4644E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2116.49786378175        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  9.8265E-01  1.1486E+00  8.0087E-01  9.2635E-01  9.4459E-01  1.0376E+00  1.1004E+00  1.1676E+00  9.1408E-01  8.5663E-01
             9.6555E-01
 PARAMETER:  8.2500E-02  2.3851E-01 -1.2205E-01  2.3494E-02  4.2999E-02  1.3693E-01  1.9571E-01  2.5495E-01  1.0167E-02 -5.4750E-02
             6.4947E-02
 GRADIENT:  -6.0121E+00  2.1452E+01 -5.5475E+00  2.1074E+01 -1.3965E+01  3.4696E+00  4.0521E+00  5.7379E+00 -7.3383E+00  4.0635E+00
             3.2337E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2117.91679257046        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  9.8618E-01  1.2294E+00  7.1304E-01  8.5620E-01  9.5357E-01  1.0287E+00  9.8130E-01  8.4542E-01  1.0417E+00  8.7973E-01
             9.6051E-01
 PARAMETER:  8.6085E-02  3.0657E-01 -2.3822E-01 -5.5256E-02  5.2452E-02  1.2825E-01  8.1123E-02 -6.7920E-02  1.4084E-01 -2.8144E-02
             5.9704E-02
 GRADIENT:  -4.1820E-02  9.8851E-01  1.5160E-01  4.2347E+00 -4.5770E+00 -2.8795E-01 -6.3536E-01  5.3508E-01  7.5482E-01  6.5053E-01
            -1.9644E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2118.41044065405        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  9.8763E-01  1.5039E+00  5.4192E-01  6.7919E-01  1.0236E+00  1.0319E+00  8.6654E-01  4.5539E-01  1.2269E+00  9.3841E-01
             9.5935E-01
 PARAMETER:  8.7549E-02  5.0809E-01 -5.1264E-01 -2.8685E-01  1.2329E-01  1.3141E-01 -4.3252E-02 -6.8660E-01  3.0451E-01  3.6428E-02
             5.8498E-02
 GRADIENT:  -4.3612E-02  5.5467E+00 -1.7661E+00  5.4094E+00 -7.2679E-01  1.2747E-01  1.3977E-01  4.8669E-01 -3.0244E-01  5.2621E-01
             2.8940E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2118.50737497319        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      885
 NPARAMETR:  9.8900E-01  1.5408E+00  5.2353E-01  6.4657E-01  1.0407E+00  1.0324E+00  8.5073E-01  2.2813E-01  1.2669E+00  9.4546E-01
             9.5846E-01
 PARAMETER:  8.8937E-02  5.3232E-01 -5.4716E-01 -3.3608E-01  1.3994E-01  1.3185E-01 -6.1662E-02 -1.3778E+00  3.3657E-01  4.3918E-02
             5.7570E-02
 GRADIENT:   2.7742E+00 -1.8990E+00  5.2428E+00 -6.2090E+00 -2.3553E+00  3.1697E-01 -1.6022E+00  2.9556E-02 -1.5456E+00 -1.6935E+00
            -3.7803E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2118.59426017827        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1061
 NPARAMETR:  9.8763E-01  1.5485E+00  5.0654E-01  6.4401E-01  1.0354E+00  1.0316E+00  8.5407E-01  1.0387E-01  1.2775E+00  9.4782E-01
             9.5853E-01
 PARAMETER:  8.7556E-02  5.3729E-01 -5.8016E-01 -3.4004E-01  1.3478E-01  1.3107E-01 -5.7737E-02 -2.1646E+00  3.4491E-01  4.6410E-02
             5.7645E-02
 GRADIENT:  -3.5059E-01 -2.1318E+00 -1.3195E+00  5.6206E-01  1.1336E-01 -1.0181E-01  4.3888E-01  3.0444E-02  2.5697E-01  3.9920E-01
             8.4085E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2118.60436675988        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1240             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8898E-01  1.5524E+00  5.0439E-01  6.4002E-01  1.0381E+00  1.0324E+00  8.5087E-01  2.4864E-02  1.2829E+00  9.4747E-01
             9.5856E-01
 PARAMETER:  8.8914E-02  5.3982E-01 -5.8440E-01 -3.4625E-01  1.3741E-01  1.3193E-01 -6.1493E-02 -3.5943E+00  3.4914E-01  4.6035E-02
             5.7672E-02
 GRADIENT:   4.7516E+02  5.3492E+02  6.8579E+00  1.3403E+02  1.2599E+01  6.7663E+01  6.7841E+00  4.3931E-03  2.0799E+01  7.0775E-01
             1.1128E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2118.60675962738        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1410
 NPARAMETR:  9.8867E-01  1.5525E+00  5.0517E-01  6.4032E-01  1.0379E+00  1.0324E+00  8.5103E-01  1.0000E-02  1.2826E+00  9.4778E-01
             9.5846E-01
 PARAMETER:  8.8602E-02  5.3989E-01 -5.8286E-01 -3.4579E-01  1.3720E-01  1.3190E-01 -6.1306E-02 -4.6151E+00  3.4892E-01  4.6363E-02
             5.7576E-02
 GRADIENT:   1.8215E+00 -3.3243E+00 -4.7880E-01 -8.1539E-01  2.1894E-01  2.3451E-01 -3.2093E-02  0.0000E+00  1.1702E-01 -3.7731E-02
             2.0384E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1410
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.2934E-04 -2.4841E-02 -2.9460E-04  2.1729E-02 -3.0350E-02
 SE:             2.9900E-02  2.3889E-02  1.3521E-04  2.4175E-02  2.2281E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9121E-01  2.9841E-01  2.9343E-02  3.6874E-01  1.7314E-01

 ETASHRINKSD(%)  1.0000E-10  1.9969E+01  9.9547E+01  1.9012E+01  2.5357E+01
 ETASHRINKVR(%)  1.0000E-10  3.5950E+01  9.9998E+01  3.4409E+01  4.4285E+01
 EBVSHRINKSD(%)  2.9125E-01  1.9796E+01  9.9606E+01  1.9930E+01  2.4048E+01
 EBVSHRINKVR(%)  5.8166E-01  3.5673E+01  9.9998E+01  3.5887E+01  4.2312E+01
 RELATIVEINF(%)  9.9362E+01  6.0446E+00  3.2066E-04  6.3641E+00  1.1317E+01
 EPSSHRINKSD(%)  3.4206E+01
 EPSSHRINKVR(%)  5.6712E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2118.6067596273811     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1199.6682264227084     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.67
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.07
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2118.607       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.89E-01  1.55E+00  5.05E-01  6.40E-01  1.04E+00  1.03E+00  8.51E-01  1.00E-02  1.28E+00  9.48E-01  9.58E-01
 


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
+        1.06E+03
 
 TH 2
+       -3.93E+00  4.29E+02
 
 TH 3
+        7.50E+00  2.50E+02  8.03E+02
 
 TH 4
+       -4.45E+00  3.11E+02 -4.31E+02  1.07E+03
 
 TH 5
+        1.53E+00 -2.39E+02 -4.66E+02  3.41E+02  6.65E+02
 
 TH 6
+        4.45E-01 -7.03E-01  3.60E+00 -2.02E+00  3.13E-01  1.85E+02
 
 TH 7
+        6.98E-01  1.08E+01 -4.71E+01  2.65E+00 -5.62E+00  1.51E-01  1.22E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.90E-01 -1.97E+01 -3.00E+01  5.45E+01 -8.46E-01 -3.50E-01  1.69E+01  0.00E+00  5.54E+01
 
 TH10
+        1.24E+00 -1.68E+01 -4.84E+01 -6.62E+00 -6.50E+01  4.69E-01  2.04E+01  0.00E+00  9.74E+00  7.63E+01
 
 TH11
+       -7.20E+00 -1.41E+01 -1.16E+01 -9.76E+00  3.33E-01  1.25E+00  5.84E+00  0.00E+00  8.06E+00  1.67E+01  4.36E+02
 
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
 #CPUT: Total CPU Time in Seconds,       28.798
Stop Time:
Thu Sep 30 01:52:24 CDT 2021
