Sat Oct 23 15:36:42 CDT 2021
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
$DATA ../../../../data/SD1/SL3/dat83.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      978
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      878
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
 RAW OUTPUT FILE (FILE): m83.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1450.65314372561        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4104E+02  1.1251E+02  8.2952E+01  1.3726E+02  4.5822E+01  3.5580E+01 -7.4149E+01 -1.4622E+02 -6.0288E+01  1.3635E+01
            -4.3420E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2819.44173322020        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0783E+00  1.1502E+00  1.0358E+00  9.2930E-01  1.1158E+00  1.0294E+00  1.0500E+00  9.6940E-01  9.2618E-01  8.9161E-01
             2.0508E+00
 PARAMETER:  1.7542E-01  2.3991E-01  1.3517E-01  2.6673E-02  2.0960E-01  1.2894E-01  1.4877E-01  6.8922E-02  2.3312E-02 -1.4722E-02
             8.1822E-01
 GRADIENT:   3.1454E+02  6.4080E+01 -3.5583E+00  2.8262E+01  1.2328E+01  1.0941E+01  2.6976E+00 -6.0453E+00 -1.8968E+01 -1.5525E+01
            -3.0486E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2828.27275092511        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0106E+00  9.5705E-01  7.5314E-01  1.0220E+00  8.8013E-01  9.9041E-01  1.2875E+00  3.4239E-01  9.6527E-01  6.8911E-01
             2.1771E+00
 PARAMETER:  1.1055E-01  5.6102E-02 -1.8350E-01  1.2177E-01 -2.7690E-02  9.0359E-02  3.5270E-01 -9.7180E-01  6.4653E-02 -2.7235E-01
             8.7801E-01
 GRADIENT:   8.7777E+01  3.2432E+01 -6.0828E+01  4.3471E+01  5.4492E+01 -1.3418E-01  1.1225E+01 -1.4556E+00 -8.5538E+00  3.7352E-01
            -1.6325E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2839.18425880246        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  9.8125E-01  1.0284E+00  9.9124E-01  9.7972E-01  1.0252E+00  1.0062E+00  1.0768E+00  9.3750E-01  9.4668E-01  8.6904E-01
             2.3260E+00
 PARAMETER:  8.1067E-02  1.2802E-01  9.1197E-02  7.9515E-02  1.2486E-01  1.0616E-01  1.7395E-01  3.5460E-02  4.5205E-02 -4.0366E-02
             9.4416E-01
 GRADIENT:   1.5133E+00 -9.0569E+00 -5.6493E+00  3.2129E+00  7.4178E+00  4.7160E+00 -1.2412E+00  1.5695E+00 -1.7445E+00 -2.3587E+00
             1.6150E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2850.73691853005        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0165E+00  1.4209E+00  1.3445E+00  8.1146E-01  1.4385E+00  1.0264E+00  8.3310E-01  1.8281E+00  1.0205E+00  1.2983E+00
             2.3082E+00
 PARAMETER:  1.1641E-01  4.5127E-01  3.9599E-01 -1.0893E-01  4.6358E-01  1.2601E-01 -8.2597E-02  7.0330E-01  1.2030E-01  3.6106E-01
             9.3646E-01
 GRADIENT:  -9.5595E+00  2.8037E+01 -1.1639E+01  6.2113E+01  3.6704E+00  5.1557E+00 -5.8658E-01 -1.7051E+00  1.0623E-01  1.3554E+01
             5.5117E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2860.13487124914        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      546
 NPARAMETR:  1.0137E+00  1.7283E+00  1.8332E+00  5.8647E-01  1.7854E+00  9.9880E-01  7.6767E-01  3.6217E+00  1.0852E+00  1.4302E+00
             2.2755E+00
 PARAMETER:  1.1363E-01  6.4715E-01  7.0605E-01 -4.3364E-01  6.7963E-01  9.8794E-02 -1.6439E-01  1.3869E+00  1.8180E-01  4.5785E-01
             9.2220E-01
 GRADIENT:  -1.5189E+01 -1.8925E+00 -6.6584E+00  2.0093E+00 -1.1150E+00 -5.3655E+00 -1.8992E+00  1.4116E+00  1.0335E+00  4.8158E-01
            -5.5680E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2860.70834888147        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      723
 NPARAMETR:  1.0271E+00  1.7271E+00  2.0736E+00  5.8729E-01  1.8152E+00  1.0079E+00  8.0432E-01  3.8515E+00  9.8417E-01  1.4608E+00
             2.2837E+00
 PARAMETER:  1.2675E-01  6.4647E-01  8.2930E-01 -4.3223E-01  6.9617E-01  1.0785E-01 -1.1775E-01  1.4485E+00  8.4040E-02  4.7901E-01
             9.2579E-01
 GRADIENT:   1.3135E+01 -1.6392E+00 -2.7087E+00 -3.9349E+00 -1.8841E+00 -1.7512E+00  1.4537E+00  6.4406E-01 -8.4814E-01  2.5019E+00
             2.8594E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2862.93001473447        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      900
 NPARAMETR:  1.0123E+00  1.3518E+00  3.3995E+00  8.5760E-01  1.6996E+00  1.0192E+00  9.4509E-01  3.6030E+00  7.6520E-01  1.2755E+00
             2.2641E+00
 PARAMETER:  1.1219E-01  4.0141E-01  1.3236E+00 -5.3617E-02  6.3041E-01  1.1903E-01  4.3525E-02  1.3818E+00 -1.6761E-01  3.4337E-01
             9.1719E-01
 GRADIENT:  -1.7081E+01  2.2410E+01 -3.0368E+00  2.1321E+01  4.3085E+00  2.2942E+00 -3.2009E+00 -3.8490E+00 -3.3861E+00  1.1324E+00
            -1.4298E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2863.79898043414        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1076
 NPARAMETR:  1.0216E+00  1.2036E+00  4.1431E+00  9.4953E-01  1.6442E+00  1.0092E+00  1.0100E+00  3.8139E+00  7.6258E-01  1.1779E+00
             2.2729E+00
 PARAMETER:  1.2142E-01  2.8534E-01  1.5215E+00  4.8213E-02  5.9725E-01  1.0916E-01  1.0998E-01  1.4387E+00 -1.7105E-01  2.6372E-01
             9.2106E-01
 GRADIENT:   2.9017E+00  7.5396E+00  1.3818E+00  5.1904E+00 -1.9046E+00 -9.8212E-01 -5.1089E-01 -3.1496E-01  9.3785E-01 -1.0683E+00
            -8.1433E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2863.91121770794        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1238
 NPARAMETR:  1.0202E+00  1.1286E+00  4.1939E+00  9.9452E-01  1.6114E+00  1.0115E+00  1.0669E+00  3.7724E+00  7.2705E-01  1.1441E+00
             2.2719E+00
 PARAMETER:  1.1996E-01  2.2097E-01  1.5336E+00  9.4507E-02  5.7712E-01  1.1144E-01  1.6474E-01  1.4277E+00 -2.1875E-01  2.3458E-01
             9.2063E-01
 GRADIENT:  -4.8547E-02 -6.4327E-02  2.5602E-02 -9.7481E-02 -5.2994E-02 -1.0416E-02 -7.0088E-02 -1.4166E-02  7.7658E-02 -6.1688E-02
            -1.4811E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1238
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5582E-03 -1.0959E-02 -4.4146E-02  4.0254E-03 -4.2460E-02
 SE:             2.9575E-02  2.1871E-02  1.9907E-02  2.0317E-02  2.1378E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5798E-01  6.1633E-01  2.6581E-02  8.4294E-01  4.7011E-02

 ETASHRINKSD(%)  9.1920E-01  2.6729E+01  3.3308E+01  3.1936E+01  2.8383E+01
 ETASHRINKVR(%)  1.8300E+00  4.6313E+01  5.5522E+01  5.3672E+01  4.8709E+01
 EBVSHRINKSD(%)  1.1316E+00  2.7078E+01  3.5523E+01  3.4229E+01  2.5379E+01
 EBVSHRINKVR(%)  2.2505E+00  4.6824E+01  5.8428E+01  5.6742E+01  4.4317E+01
 RELATIVEINF(%)  9.7701E+01  6.1082E+00  2.1392E+01  5.0167E+00  3.0925E+01
 EPSSHRINKSD(%)  1.7988E+01
 EPSSHRINKVR(%)  3.2740E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          878
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1613.6560643074051     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2863.9112177079378     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1250.2551534005327     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.64
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2863.911       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.13E+00  4.19E+00  9.95E-01  1.61E+00  1.01E+00  1.07E+00  3.77E+00  7.27E-01  1.14E+00  2.27E+00
 


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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       78.731
Stop Time:
Sat Oct 23 15:36:55 CDT 2021
