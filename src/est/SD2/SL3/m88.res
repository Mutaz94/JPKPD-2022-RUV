Sat Oct 23 19:32:38 CDT 2021
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
$DATA ../../../../data/SD2/SL3/dat88.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      799
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

 TOT. NO. OF OBS RECS:      699
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
 RAW OUTPUT FILE (FILE): m88.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1790.00266439253        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1837E+02 -1.2534E+02 -3.7134E+01  1.3313E+02  1.8323E+02  4.9478E+01 -2.7647E+01 -6.6310E+01 -4.9628E+01  3.3862E+01
            -2.2592E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2444.65055243774        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.7782E-01  1.2733E+00  1.0176E+00  8.5219E-01  1.0250E+00  8.6315E-01  8.7939E-01  1.0540E+00  1.1162E+00  6.2647E-01
             1.9436E+00
 PARAMETER:  7.7572E-02  3.4164E-01  1.1743E-01 -5.9946E-02  1.2468E-01 -4.7163E-02 -2.8529E-02  1.5258E-01  2.0991E-01 -3.6765E-01
             7.6452E-01
 GRADIENT:   4.0848E+00  3.9998E+01  1.1737E+01 -1.3656E+01  2.2690E+01 -4.5655E+01 -6.9665E+00 -1.2743E+01 -1.8168E+01 -1.6894E+01
            -5.3986E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2455.90987882711        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      245
 NPARAMETR:  1.0058E+00  1.7343E+00  6.9949E-01  6.0156E-01  1.1739E+00  9.3478E-01  8.1118E-01  1.4054E+00  1.6500E+00  8.4879E-01
             1.9279E+00
 PARAMETER:  1.0577E-01  6.5061E-01 -2.5740E-01 -4.0823E-01  2.6032E-01  3.2557E-02 -1.0926E-01  4.4029E-01  6.0078E-01 -6.3938E-02
             7.5643E-01
 GRADIENT:  -2.6194E+01  3.0565E+01 -1.5671E+01  2.2830E+01 -1.6406E+01 -1.9864E+01  2.0441E+01  4.5938E+00  4.3503E+00 -6.0621E+00
            -5.6962E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2462.17733102200        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      421
 NPARAMETR:  1.0116E+00  1.7120E+00  9.8754E-01  6.1675E-01  1.2758E+00  9.7148E-01  6.7920E-01  2.1148E+00  1.5815E+00  9.4860E-01
             1.9930E+00
 PARAMETER:  1.1149E-01  6.3768E-01  8.7465E-02 -3.8329E-01  3.4355E-01  7.1060E-02 -2.8683E-01  8.4894E-01  5.5837E-01  4.7227E-02
             7.8966E-01
 GRADIENT:  -1.2683E+01  4.7530E+00 -9.5473E+00  1.6979E+01  1.8957E+01 -3.3916E+00  4.2466E+00  1.3027E+00 -4.4166E-01 -1.3734E+00
             6.5997E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2465.91671049096        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      596
 NPARAMETR:  1.0190E+00  1.5326E+00  1.9225E+00  7.1767E-01  1.3020E+00  9.7124E-01  5.1127E-01  3.0924E+00  1.4615E+00  9.7006E-01
             1.9824E+00
 PARAMETER:  1.1878E-01  5.2694E-01  7.5365E-01 -2.3174E-01  3.6391E-01  7.0821E-02 -5.7086E-01  1.2289E+00  4.7945E-01  6.9605E-02
             7.8433E-01
 GRADIENT:   6.6090E+00  9.8592E-01  4.1100E+00 -1.3547E+01 -1.8007E+01 -2.5419E+00  2.4504E+00 -2.9631E-01 -3.5595E+00  4.1826E+00
             3.0096E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2467.47593321138        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      771
 NPARAMETR:  1.0193E+00  1.2549E+00  2.3249E+00  9.2355E-01  1.2436E+00  9.7830E-01  4.0742E-01  3.1015E+00  1.2655E+00  8.3828E-01
             1.9798E+00
 PARAMETER:  1.1907E-01  3.2702E-01  9.4368E-01  2.0468E-02  3.1803E-01  7.8058E-02 -7.9790E-01  1.2319E+00  3.3547E-01 -7.6399E-02
             7.8298E-01
 GRADIENT:   7.7619E+00  6.2337E+00  4.2919E-01  9.7407E+00  3.2394E+00  1.7402E-01  2.0893E+00  7.3635E-01  5.4742E-02  1.5376E-01
            -3.5515E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2468.45311786995        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      947
 NPARAMETR:  1.0172E+00  1.2495E+00  2.2677E+00  9.0762E-01  1.2370E+00  9.8051E-01  1.6349E-01  3.1233E+00  1.3270E+00  8.3689E-01
             1.9767E+00
 PARAMETER:  1.1707E-01  3.2273E-01  9.1876E-01  3.0672E-03  3.1267E-01  8.0319E-02 -1.7110E+00  1.2389E+00  3.8289E-01 -7.8065E-02
             7.8144E-01
 GRADIENT:   3.7577E+00 -4.6081E+00 -8.4107E-01 -3.8967E+00  2.2081E+00  1.2340E+00  3.7364E-01  2.7329E-01 -7.0168E-02 -1.1960E+00
            -2.4112E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2468.64511073994        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1123
 NPARAMETR:  1.0158E+00  1.2409E+00  2.3191E+00  9.1527E-01  1.2392E+00  9.7740E-01  3.4462E-02  3.1330E+00  1.3323E+00  8.4695E-01
             1.9801E+00
 PARAMETER:  1.1565E-01  3.1583E-01  9.4116E-01  1.1464E-02  3.1450E-01  7.7142E-02 -3.2679E+00  1.2420E+00  3.8689E-01 -6.6109E-02
             7.8316E-01
 GRADIENT:   4.0782E-01 -6.1106E-01 -1.8844E-01  3.9845E-02  8.1252E-01  3.8959E-02  1.8033E-02  4.1466E-03  6.9578E-02 -9.0703E-02
             2.5980E-01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -2468.65255042126        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1250
 NPARAMETR:  1.0156E+00  1.2435E+00  2.3235E+00  9.1353E-01  1.2402E+00  9.7719E-01  1.0000E-02  3.1382E+00  1.3351E+00  8.4940E-01
             1.9800E+00
 PARAMETER:  1.1544E-01  3.1795E-01  9.4305E-01  9.5601E-03  3.1531E-01  7.6923E-02 -4.6562E+00  1.2436E+00  3.8904E-01 -6.3226E-02
             7.8312E-01
 GRADIENT:  -1.0198E-01  2.4319E-01  2.4885E-02  2.9504E-01  8.9931E-02 -4.7282E-02  0.0000E+00 -7.7977E-03 -4.8056E-02  2.8763E-02
             2.3257E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1250
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.8903E-03 -1.1483E-03 -3.8465E-02  2.1969E-03 -3.1435E-02
 SE:             2.9689E-02  3.0425E-04  1.9452E-02  2.8680E-02  2.0454E-02
 N:                     100         100         100         100         100

 P VAL.:         9.2245E-01  1.6059E-04  4.7992E-02  9.3894E-01  1.2432E-01

 ETASHRINKSD(%)  5.3776E-01  9.8981E+01  3.4833E+01  3.9198E+00  3.1478E+01
 ETASHRINKVR(%)  1.0726E+00  9.9990E+01  5.7533E+01  7.6859E+00  5.3048E+01
 EBVSHRINKSD(%)  1.0437E+00  9.9143E+01  3.6467E+01  4.0048E+00  3.0296E+01
 EBVSHRINKVR(%)  2.0766E+00  9.9993E+01  5.9635E+01  7.8492E+00  5.1413E+01
 RELATIVEINF(%)  9.7868E+01  1.2960E-03  1.8538E+01  1.9581E+01  2.1258E+01
 EPSSHRINKSD(%)  2.1900E+01
 EPSSHRINKVR(%)  3.9004E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          699
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1284.6760694201323     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2468.6525504212632     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1183.9764810011309     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.67
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2468.653       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.24E+00  2.32E+00  9.14E-01  1.24E+00  9.77E-01  1.00E-02  3.14E+00  1.34E+00  8.49E-01  1.98E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.01
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,      111.511
Stop Time:
Sat Oct 23 19:32:55 CDT 2021
