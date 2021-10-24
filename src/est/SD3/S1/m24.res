Sat Oct 23 22:57:13 CDT 2021
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
$DATA ../../../../data/SD3/S1/dat24.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m24.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2116.94276421490        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5804E+02 -3.4317E+01  1.7458E+01 -6.1166E+01 -2.1010E+01  3.3198E+01 -3.7213E+00  6.7408E+00  9.2957E-01 -3.0621E-01
             4.1088E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2124.83109816431        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.7275E-01  1.0415E+00  1.0071E+00  1.0609E+00  1.0235E+00  1.0130E+00  1.0236E+00  9.7394E-01  9.8696E-01  1.0158E+00
             9.3875E-01
 PARAMETER:  7.2375E-02  1.4063E-01  1.0707E-01  1.5907E-01  1.2318E-01  1.1293E-01  1.2335E-01  7.3598E-02  8.6875E-02  1.1572E-01
             3.6789E-02
 GRADIENT:   1.1589E+00  4.8112E+00  5.5735E+00 -5.0783E-01 -1.3161E+01  4.2092E-01 -6.0965E+00  4.2022E+00  2.3033E+00 -3.4541E+00
            -1.0015E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2126.40016715944        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  9.7212E-01  9.2263E-01  1.0066E+00  1.1458E+00  9.8856E-01  1.0417E+00  1.4583E+00  7.1993E-01  8.1560E-01  1.0336E+00
             9.4777E-01
 PARAMETER:  7.1719E-02  1.9474E-02  1.0656E-01  2.3611E-01  8.8490E-02  1.4090E-01  4.7729E-01 -2.2861E-01 -1.0383E-01  1.3300E-01
             4.6359E-02
 GRADIENT:   2.6830E+00  2.1757E+01  6.4007E+00  2.2504E+01  4.8038E+00  1.1650E+01  4.4843E+00 -4.1831E+00 -7.2308E+00  1.5003E+00
            -1.2965E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2127.88399785004        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  9.7245E-01  8.9950E-01  8.0339E-01  1.1263E+00  8.4942E-01  1.0116E+00  1.3831E+00  4.7446E-01  8.6831E-01  8.7624E-01
             9.4723E-01
 PARAMETER:  7.2065E-02 -5.9114E-03 -1.1892E-01  2.1898E-01 -6.3197E-02  1.1149E-01  4.2435E-01 -6.4558E-01 -4.1211E-02 -3.2110E-02
             4.5783E-02
 GRADIENT:  -5.1796E-02  4.5346E+00  2.1830E+00  4.3443E+00 -3.6522E+00 -1.8169E-01 -3.5957E-01 -6.0388E-01  4.2613E-01  3.4965E-01
            -7.2338E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2128.33889180955        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      709
 NPARAMETR:  9.6775E-01  6.4844E-01  1.1074E+00  1.3054E+00  9.1629E-01  1.0064E+00  1.7959E+00  8.9146E-01  7.6674E-01  9.6449E-01
             9.4999E-01
 PARAMETER:  6.7222E-02 -3.3319E-01  2.0203E-01  3.6648E-01  1.2577E-02  1.0636E-01  6.8549E-01 -1.4896E-02 -1.6561E-01  6.3844E-02
             4.8699E-02
 GRADIENT:   3.4064E-02  8.6605E+00  4.9578E+00  1.2059E+01 -8.9488E+00 -3.9844E-01 -7.5161E-01  4.7368E-01 -1.7050E+00 -2.3896E-01
             5.7599E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2128.86547209907        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      887
 NPARAMETR:  9.6293E-01  3.9316E-01  1.4820E+00  1.4818E+00  1.0055E+00  1.0033E+00  2.5005E+00  1.2367E+00  7.0794E-01  1.0794E+00
             9.4900E-01
 PARAMETER:  6.2226E-02 -8.3354E-01  4.9337E-01  4.9326E-01  1.0547E-01  1.0333E-01  1.0165E+00  3.1243E-01 -2.4540E-01  1.7636E-01
             4.7649E-02
 GRADIENT:  -4.6801E-01  6.4618E+00  3.1711E+00  1.6606E+01 -2.9075E+00 -9.2208E-03  4.9987E-01  7.5484E-02 -1.8588E+00 -3.3808E-01
            -1.2369E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2129.13434482289        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1069
 NPARAMETR:  9.6245E-01  3.1852E-01  1.5524E+00  1.5193E+00  1.0152E+00  1.0025E+00  2.7702E+00  1.3033E+00  7.0367E-01  1.0963E+00
             9.5033E-01
 PARAMETER:  6.1724E-02 -1.0441E+00  5.3983E-01  5.1827E-01  1.1510E-01  1.0245E-01  1.1189E+00  3.6490E-01 -2.5144E-01  1.9198E-01
             4.9051E-02
 GRADIENT:   1.1035E+00  1.1559E+00  1.0362E+00 -6.8227E+00  1.7160E+00  1.2811E-01  4.9004E-02  3.7921E-01 -1.8161E-01  1.3035E-01
             1.0082E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2129.17117253493        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1244
 NPARAMETR:  9.6102E-01  2.8053E-01  1.5281E+00  1.5468E+00  9.9114E-01  1.0009E+00  2.8999E+00  1.2786E+00  7.0670E-01  1.0785E+00
             9.5015E-01
 PARAMETER:  6.0238E-02 -1.1711E+00  5.2401E-01  5.3621E-01  9.1101E-02  1.0093E-01  1.1647E+00  3.4573E-01 -2.4715E-01  1.7557E-01
             4.8867E-02
 GRADIENT:  -7.1503E-01  1.7288E+00  9.1671E-01  5.6712E+00 -1.9251E-01 -1.6320E-01 -6.0699E-01 -6.0084E-01 -5.0108E-01  3.3794E-01
            -2.9505E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2129.17401160415        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1426
 NPARAMETR:  9.6240E-01  2.6384E-01  1.5216E+00  1.5391E+00  9.8718E-01  1.0016E+00  3.0925E+00  1.2891E+00  7.1083E-01  1.0710E+00
             9.5052E-01
 PARAMETER:  6.1677E-02 -1.2324E+00  5.1975E-01  5.3119E-01  8.7096E-02  1.0164E-01  1.2290E+00  3.5394E-01 -2.4132E-01  1.6855E-01
             4.9256E-02
 GRADIENT:   3.1654E+00 -1.5230E+00 -3.0389E-01 -3.2140E+01  3.0219E+00  2.0812E-01  2.3586E+00  9.5531E-01  4.3617E+00  1.3119E-02
             6.5404E-01

0ITERATION NO.:   44    OBJECTIVE VALUE:  -2129.26897208165        NO. OF FUNC. EVALS.: 132
 CUMULATIVE NO. OF FUNC. EVALS.:     1558
 NPARAMETR:  9.6181E-01  2.6842E-01  1.5219E+00  1.5455E+00  9.8517E-01  1.0014E+00  3.0536E+00  1.2812E+00  7.0079E-01  1.0709E+00
             9.5024E-01
 PARAMETER:  6.1061E-02 -1.2152E+00  5.1997E-01  5.3535E-01  8.5062E-02  1.0144E-01  1.2163E+00  3.4779E-01 -2.5554E-01  1.6848E-01
             4.8959E-02
 GRADIENT:  -2.3288E-02  1.4657E-01  1.1736E+00  1.9218E+00 -2.7062E-01 -3.2583E-02 -2.1535E-01  1.2770E-01  2.1043E-01  8.0612E-02
             1.0182E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1558
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3107E-03  3.2117E-02 -4.2512E-02 -2.8354E-02 -3.4405E-02
 SE:             2.9899E-02  1.7107E-02  1.8245E-02  2.4591E-02  2.0302E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6503E-01  6.0468E-02  1.9803E-02  2.4889E-01  9.0148E-02

 ETASHRINKSD(%)  1.0000E-10  4.2688E+01  3.8877E+01  1.7618E+01  3.1985E+01
 ETASHRINKVR(%)  1.0000E-10  6.7153E+01  6.2640E+01  3.2133E+01  5.3739E+01
 EBVSHRINKSD(%)  3.4368E-01  5.0987E+01  4.1916E+01  1.2742E+01  2.7770E+01
 EBVSHRINKVR(%)  6.8619E-01  7.5977E+01  6.6263E+01  2.3860E+01  4.7828E+01
 RELATIVEINF(%)  9.8739E+01  4.3648E+00  9.7706E+00  1.4569E+01  1.5183E+01
 EPSSHRINKSD(%)  3.4532E+01
 EPSSHRINKVR(%)  5.7140E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2129.2689720816534     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1210.3304388769807     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.22
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2129.269       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.62E-01  2.68E-01  1.52E+00  1.55E+00  9.85E-01  1.00E+00  3.05E+00  1.28E+00  7.01E-01  1.07E+00  9.50E-01
 


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
 #CPUT: Total CPU Time in Seconds,      126.700
Stop Time:
Sat Oct 23 22:57:32 CDT 2021
