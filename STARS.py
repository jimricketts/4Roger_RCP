# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 16:06:31 2014

@author: s4493222
"""
'''
'Sequential regime shift detection in the mean
'with red noise removal
' Version 3.2
'Author: Sergei Rodionov
'Date modified: 04/26/2006
'''
import numpy as np
import scipy.stats as stats
import os
SVNRevision="$Revision: 286 $"

fn=os.environ["HOMEPATH"]+"\\Documents\\abrupt\\CMIP5_breakpoints\\Artyfishul_dartar4Jim (1).csv"

alpha = np.NaN #crashout unless properly initialised

def EqN(Ns, alpha):
  '''
    'calculates the equivalent sample size
    'as in von Storch and Zwiers (1999, p. 115)
    'Note: The formula for the equiv. sample size published by
    'Zwiers and Storch (1995)in J. Climate is somewhat different
    'Ns original sample size
    'alpha - AR1 coefficient defined for the entire module
    sum = sum + (1 - i / Ns) * alpha ^ i
  '''
  suma =0.0
  for i in range(1, Ns):
    suma += (1.0 - i / float(Ns)) * alpha ** i
  return max(2, min(Ns, Ns / (1 + suma)))
  
def EqP(rng1, rng2, alpha):
  '''
    calculates t-test using the equivalent sample size
    'as in von Storch and Zwiers (1999, p. 115)
    'Note: The formula for the equiv. sample size published in Zwiers and Storch (1995)
    'in J. Climate is somewhat different
    #two ranges of values are compared - since python works with slices, rather than computing slices internally as in VBA code
  '''
  #calculate the means and variances
  ave1= np.mean(rng1, axis=0)
  ave2 =np.mean(rng2, axis =0) 
  var1= np.var(rng1, axis=0)
  var2 =np.var(rng2, axis =0) 
  #'calculate effective sample sizes
  #'range1
  Ns = len(rng1)
  if Ns < 2:
    raise Exception("There is no data to calculate the effective sample size Range 1 not specified")
  eN1 = EqN(Ns, alpha)
  Ns = len(rng2)
  if Ns < 2:
    raise Exception("There is no data to calculate the effective sample size Range 2 not specified")
  eN2 = EqN(Ns, alpha)
   # t-statistics
  T_stat = np.sqrt(var1 / eN1 + var2 / eN2)
  T_stat = np.abs(ave1 - ave2) / T_stat
   #probability for the t-distribution 2 tailed
  #return T_stat,  eN1 + eN2 - 2, 2.0 * (1.0-stats.t.cdf(T_stat,  eN1 + eN2 - 2))
  return 2.0 * (1.0-stats.t.cdf(T_stat,  eN1 + eN2 - 2))
  #return  WorksheetFunction.TDist(T_stat, eN1 + eN2 - 2, 2)

def MPK(est, Nsub):
  '''
    'Calculates bias-corrected AR(1) coefficient using the formula of
    'Marriott-Pope and Kendall(see Orcutt and Winokur, 1969, Econometrica, 37:1,1-14)
    'est is the OLS estimate of AR(1)
    'Nsub - sample size declared globaly
  '''
  if Nsub > 4:
    result = ((float(Nsub) - 1.0) * float(est) + 1.0) / (Nsub - 4.0)
  else:
    #'should not be here
    result = est
  return result

def IPN4(est, Nsub):
  '''
    'Calculates bias-corrected AR(1) coefficient [est] using
    'interative 4-step procedure with corrections that are
    'reverse proportional to the subsample size - nsub
  '''
  result = est + 1 / float(Nsub)
  for i in range(3):
    result += np.abs(result) / float(Nsub)
  return result

def OLSAR1(x):
  '''
    'Function OLSAR1(x() As Double) As Double
    'OLS estimate of AR1 coefficient
    'from Orcutt and Winokur, 1969, Econometrica, 37:1,1-14
    Ordinary Least Squares Estimate of AR1
  '''
  ave1=np.mean(x[1:])
  ave2=np.mean(x[:-1])
  part1 = x[1:] - ave1
  part2 = x[:-1] - ave2
  sumNom=np.sum(part1 * part2)
  sumDenom=np.sum(part2 * part2)
  if sumDenom > 0.0:
    return sumNom / sumDenom
  else:
    return 0.0

def AlphaEst(x, Nsub, option=None, returnmsgs=True): # As Double
  '''
    'Public Function AlphaEst(x() As Double) As Double
    'Calculate the bias-corrected estimate of alpha
    'in x(t)= alpha * x(t-1) + white noise as the following:
    ' 1. Calculate est1, the OLS estimate of alpha for all
    '   sequential subsamples of size nsub;
    ' 2. Let est1 be the median value of all est1 for subsamples
    ' 3. Make bias correction for est1
  '''
  msg = []
  if option != None and not option in ["optRNKendall", "optIPN4"]:
    raise Exception("option nust be one of [None, 'optRNKendall', 'optIPN4']")
  elif option == "optRNKendall" and  Nsub < 5:
    Nsub = 5
    msg.append("The subsample size is too small for optRNKendall and changed to 5.")
  #'test the subsample size (Nsub) limits
  elif option == "optIPN4" and Nsub < 3:
    Nsub = 3
    msg.append("The subsample size is too small for optIPN4 and changed to 3.")
  N=len(x)
  
  if Nsub > N:
    Nsub = N
    msg.append("The subsample size greater than N, set to N ("+str(N)+")")
  #'Subsampling
  ss=np.array([OLSAR1(x[i:i+Nsub + 1]) for i in range(N - Nsub + 1)])
  est1 = np.median(ss)
  #print ss
#' bias correction
  if option == None:
    result= est1
  elif option == "optRNKendall":
    result=MPK(est1, Nsub)
  else:
    result=IPN4(est1, Nsub)
#'limits
  if returnmsgs:
    return min(0.99, max(0.0, result)), msg
  else:
    return min(0.99, max(0.0, result))

def WeighedAverage(dblW, HuberParam):
  '''
    #Function WeighedAverage(dblW() As Double) As Double
    'Calculation of the mean estimate for a given range using Huber's weights
    'The Huber's function is defined as min( 1, param/(|res|/scale))
    'param: a given parameter (txtHuberParam); affects the range where weights = 1
    'res: deviation from the mean estimate (I used the median as the first approximation)
    'scale: estimate of variation. Some use 1.483*(median absolute deviation, or MAD,
    ' of the deviations of the data from their median). I used std, i.e. dblAveVar
    'See http://www.stat.berkeley.edu/~stark/Preprints/Oersted/writeup.htm
    'The procedure consists of two iterations. At the first iteration the estimate of the
    'regime average is the simple unweighed arithmetic mean. In the second iteration
    'it is weighed average from the first iteration.
  '''
  dblEstAve = np.mean(dblW, axis=0)
  dblAveStd=np.std(dblW, axis=0)
  for k in range(2):
    dblXdev = (dblW - dblEstAve) / dblAveStd
    dblXweight = np.where(dblXdev == 0.0, 1.0, np.array([min(1,x) for  x in (HuberParam / np.abs(dblXdev))]))
    dblSumofWeights=np.sum(dblXweight,axis=0)
    dblSumAve=np.sum(dblXdev * dblXweight,axis=0)
    dblSumAve /= dblSumofWeights
    dblSumAve = dblSumAve * dblAveStd + dblEstAve
    dblEstAve = dblSumAve
  return dblEstAve

#def cusumUp(HuberParam):
#  '''
#    'This function calculates cusum for upward shifts
#  '''
#Dim i As Integer
#Dim dblXdev As Double 'normalized deviation of x value from the expected value of the new
#                        'regime, i.e., from dblRegimeMean + dblDiff
#Dim dblXweight As Double 'Huber's weight of value x.
#Dim dblSumofWeights As Double
#                        
#cusumUp = 0
#For i = intYear To intYear + cutoff - 1
#    If i > N Then
#        If dblSumofWeights > 0 Then
#            cusumUp = cusumUp / dblSumofWeights
#            Exit Function
#        End If
#    End If
#    
#    dblXdev = (dblX(i) - dblRegimeMean - dblDiff) / Sqr(dblAveVar)
#    
#    'determine the weight of the normalized deviation
#    If dblXdev = 0 Then
#        dblXweight = 1
#    Else
#        dblXweight = WorksheetFunction.Min(1, frmMain.txtHuberParam.Value / Abs(dblXdev))
#    End If
#    
#    'sum weights and weighed values
#    dblSumofWeights = dblSumofWeights + dblXweight
#    cusumUp = cusumUp + dblXdev * dblXweight
#    
#    'check if cusum turns zero
#        If cusumUp < 0 Then
#        cusumUp = 0
#        Exit Function
#    End If
#Next i
#cusumUp = cusumUp / dblSumofWeights
#
#End Function
#Function cusumDown() As Double
#'This function calculates cusum fo downward shifts
#Dim i As Integer
#Dim dblXdev As Double 'normalized deviation of x value from the expected value of the new
#                        'regime, i.e., from dblRegimeMean + dblDiff
#Dim dblXweight As Double 'Huber's weight of value x.
#Dim dblSumofWeights As Double
#
#cusumDown = 0
#For i = intYear To intYear + cutoff - 1
#    If i > N Then
#        If dblSumofWeights > 0 Then
#            cusumDown = cusumDown / dblSumofWeights
#            Exit Function
#        End If
#    End If
#    
#        dblXdev = (dblX(i) - dblRegimeMean + dblDiff) / Sqr(dblAveVar)
#    
#    'determine the weight of the normalized deviation
#    If dblXdev = 0 Then
#        dblXweight = 1
#    Else
#        dblXweight = WorksheetFunction.Min(1, frmMain.txtHuberParam.Value / Abs(dblXdev))
#    End If
#    
#    'sum weights and weighed values
#    dblSumofWeights = dblSumofWeights + dblXweight
#    cusumDown = cusumDown + dblXdev * dblXweight
#    
#    'check if cusum turns zero
#    If cusumDown > 0 Then
#        cusumDown = 0
#        Exit Function
#    End If
#Next i
#cusumDown = cusumDown / dblSumofWeights
#
#End Function
'''

Option Explicit
Public alpha As Double 'AR(1) estimate
Public Nsub As Integer 'size of subsamples

'Public prob As Double 'probability level - declared in basShiftInVariance
'Public cutoff As Integer 'cutoff length - declared in basShiftInVariance

Dim dblX() As Double 'Individual time series
Dim rsi() As Double 'cumulative sum of squares
Dim N As Integer 'total number of observations in dblX and rsi
Dim intYear As Integer
Dim ChangePoint As Integer
Dim dblDiff As Double 'critical difference between 2 regimes
Dim dblRegimeMean As Double 'mean velue of the current regime
Dim dblAveVar As Double 'average variance for the cutoff period in the series

Sub ShiftMeanMain(rngMatrix As Range)

Dim intMaxRows As Integer
Dim intMaxCols As Integer 'Years are in the first column
Dim intStartRow As Integer
Dim intEndRow As Integer
Dim intTS As Integer 'counter of time series
Dim intRow As Integer 'counter of rows
Dim intRow1 As Integer 'also a counter
Dim intCol As Integer 'a counter for columns
Dim shtNew As Worksheet
Dim shtSummary As Worksheet
Dim strDataShtName As String
Dim wbkOut As Workbook ' workbook for output
Dim strNewShtName As String
Dim intRegimeStart As Integer
Dim intRegimeEnd As Integer
Dim blnEmpty As Boolean
Dim dblMean() As Double 'use to calc the mean for selected regimes
Dim dblMeanRegime As Double 'mean value for each regime
Dim dblMeanRegimeW As Double 'weighed mean for each regime
Dim sngResiduals() As Single
Dim strActiveWbkName As String 'to remember active workbook
Dim intCountShifts As Integer 'Count the number of shifts in each year in the entire set in the summary sheet
Dim i As Integer ' a counter
Dim j As Integer 'a counter
Dim rngW As Range 'working array
Dim n1 As Integer 'used to calculate the length of the regime
Dim intDT As Integer '=1 if prewhitening is used, =0 otherwise
Dim intLagF As Integer '=1 if return to original, =0 otherwise

intMaxRows = rngMatrix.Rows.Count
intMaxCols = rngMatrix.Columns.Count
If intMaxCols > 251 Then
    MsgBox "The number of columns in the data matrix exceeds the limit of 250.", vbCritical, _
    "Too Many Columns"
    Exit Sub
End If
If Not IsNumeric(rngMatrix.Cells(intMaxRows, 1)) Then
    intMaxRows = intMaxRows - 2
End If

'status of prewhitening and the type output
If frmMain.chkPrewhitening.Value = True Then
    intDT = 1
    If frmMain.chkFilteredData = True Then
        intLagF = 0
    Else
        intLagF = 1
    End If
Else
    intDT = 0
    intLagF = 0
End If

'to remember where the data are
strActiveWbkName = ActiveWorkbook.Name
strDataShtName = ActiveSheet.Name

'determine a workbook for output
If frmMain.optOuputNewWkbk.Value Then
    Set wbkOut = Workbooks.Add
Else
    Set wbkOut = Application.ActiveWorkbook
End If

' Determine if the residual worksheet exists
For Each shtNew In wbkOut.Worksheets
        If shtNew.Name = "Residuals" Then GoTo Residuals_exists
    Next shtNew
    Set shtNew = wbkOut.Worksheets.Add
    shtNew.Name = "Residuals"
Residuals_exists:

' Determine if the summary worksheet exists
    For Each shtNew In wbkOut.Worksheets
        If shtNew.Name = "Summary" Then GoTo Summary_exists
    Next shtNew
    Set shtNew = wbkOut.Worksheets.Add
    shtNew.Name = "Summary"
Summary_exists:

'cycle for each time series
For intTS = 2 To intMaxCols

    'Return back to the data worksheet
    Application.Workbooks(strActiveWbkName).Worksheets(strDataShtName).Activate
    
    ' Determine where the time series starts
    For intRow = 1 To intMaxRows
        If IsNumeric(rngMatrix(intRow, intTS)) And _
            Not IsEmpty(rngMatrix(intRow, intTS)) Then
            intStartRow = intRow
            GoTo EndRow
        End If
    Next intRow
    
EndRow:
    ' Determine where the time series ends
    For intRow = intStartRow + 1 To intMaxRows
        If IsEmpty(rngMatrix(intRow, intTS)) Then
           intEndRow = intRow - 1
           GoTo TS
        End If
    Next intRow
    intEndRow = intMaxRows

TS:
    ' Form the time series
    N = intEndRow - intStartRow + 1
    ReDim dblX(1 To N)
    For i = intStartRow To intEndRow
        dblX(i - intStartRow + 1) = rngMatrix.Cells(i, intTS)
    Next i
    
    If frmMain.optRNnone.Value = False Then
        'Estimate alpha  - AR(1) value
        alpha = AlphaEst(dblX)
    End If

    If intDT = 1 Then
        'use prewhitening to remove red noise x(t)=x(t)-alpha * x(t-1)
        N = N - 1
        ReDim dblX(1 To N)
        For i = intStartRow + 1 To intEndRow
            dblX(i - intStartRow) = rngMatrix.Cells(i, intTS) - _
            alpha * rngMatrix.Cells(i - 1, intTS)
        Next i
        Call ShiftMean
        If frmMain.chkFilteredData.Value = False Then
            'get back to the original time series
            N = N + 1
            ReDim dblX(1 To N)
            'use dblX for temporary storage of rsi
            dblX(1) = 0
            For i = 2 To N
                dblX(i) = rsi(i - 1)
            Next i
            ReDim rsi(1 To N)
            For i = 1 To N
                rsi(i) = dblX(i)
            Next i
            For i = intStartRow To intEndRow
                dblX(i - intStartRow + 1) = rngMatrix.Cells(i, intTS)
            Next i
        End If
    Else
        Call ShiftMean
    End If
            
    'output
    If frmMain.chkSummaryOnly Then GoTo Summary_Sheet
    ReDim sngResiduals(1 To N)
    'add or activate a new sheet
    strNewShtName = Application.Workbooks(strActiveWbkName).Worksheets(strDataShtName).Cells(1, intTS)
    For Each shtNew In wbkOut.Worksheets
        If shtNew.Name = strNewShtName Then GoTo Sheet_exists
    Next shtNew
    Set shtNew = wbkOut.Worksheets.Add
    shtNew.Name = strNewShtName
Sheet_exists:
    shtNew.Activate
    
    'print labels in the first row
    ActiveSheet.UsedRange.Select
    Selection.Clear
    Range("a1").Select
    Cells(1, 1) = Application.Workbooks(strActiveWbkName).Worksheets(strDataShtName).Cells(1, 1)
    Cells(1, 2) = Application.Workbooks(strActiveWbkName).Worksheets(strDataShtName).Cells(1, intTS)
    Cells(1, 3) = "RSI"
    Cells(1, 4) = "Mean"
    Cells(1, 5) = "Weighed"
    Cells(1, 6) = "Length"
    Cells(1, 7) = "Conf"
    Cells(1, 8) = "Outliers"
    
    'print time series
    intNumberOfRegimes = 1
    intRegimeStart = 1
    For intRow = 1 To N
        Cells(intRow + 1, 1) = rngMatrix(intRow + intStartRow - 1 + intDT - intLagF, 1)
        Cells(intRow + 1, 2) = dblX(intRow)
        Cells(intRow + 1, 3) = rsi(intRow)
        
        If Abs(rsi(intRow)) > 0 Then
            intNumberOfRegimes = intNumberOfRegimes + 1
            'Calculate the mean for selected regimes
            intRegimeEnd = intRow - 1
            n1 = intRegimeEnd - intRegimeStart + 1
            ReDim dblMean(1 To n1)
            For i = intRegimeStart To intRegimeEnd
                dblMean(i - intRegimeStart + 1) = dblX(i)
            Next i
            dblMeanRegime = Application.WorksheetFunction.Average(dblMean)
            dblMeanRegimeW = WeighedAverage(dblMean)
            
            'print the average and calculate residuals
            For intRow1 = intRegimeStart To intRow - 1
               sngResiduals(intRow1) = dblX(intRow1) - dblMeanRegimeW
               Cells(intRow1 + 1, 4) = dblMeanRegime
               Cells(intRow1 + 1, 5) = dblMeanRegimeW
               Cells(intRow1 + 1, 6) = intRegimeEnd - intRegimeStart + 1 'to remember regime length
               If Abs(sngResiduals(intRow1)) / Sqr(dblAveVar) - frmMain.txtHuberParam.Value > 0 Then
                    Cells(intRow1 + 1, 8) = Sqr(dblAveVar) * frmMain.txtHuberParam.Value / Abs(sngResiduals(intRow1))
                    Cells(intRow1 + 1, 8).NumberFormat = "0.00"
               End If
            Next intRow1
            intRegimeStart = intRow
        End If
    Next intRow
    
    'calculate the mean for the last regime
    n1 = N - intRegimeStart + 1
    ReDim dblMean(1 To n1)
    For i = intRegimeStart To N
        dblMean(i - intRegimeStart + 1) = dblX(i)
    Next i
    dblMeanRegime = Application.WorksheetFunction.Average(dblMean)
    dblMeanRegimeW = WeighedAverage(dblMean)
    'print the average and calculate residuals for the last regime
    For intRow1 = intRegimeStart To N
        sngResiduals(intRow1) = dblX(intRow1) - dblMeanRegimeW
        Cells(intRow1 + 1, 4) = dblMeanRegime
        Cells(intRow1 + 1, 5) = dblMeanRegimeW
        Cells(intRow1 + 1, 6) = N - intRegimeStart + 1 'to remember regime length
        If Abs(sngResiduals(intRow1)) / Sqr(dblAveVar) - frmMain.txtHuberParam.Value > 0 Then
            Cells(intRow1 + 1, 8) = Sqr(dblAveVar) * frmMain.txtHuberParam.Value / Abs(sngResiduals(intRow1))
            Cells(intRow1 + 1, 8).NumberFormat = "0.00"
        Else
            Cells(intRow1 + 1, 8) = ""
        End If
    Next intRow1
    
    'calculate the real significance level
    intRegimeStart = 1
    For intRow = 2 To N - 1
        If Abs(rsi(intRow)) > 0 Then
            If intRow - intRegimeStart > 1 Then
                If frmMain.optRNnone.Value = True Or _
                (frmMain.chkPrewhitening.Value = True And frmMain.chkFilteredData.Value = True) Then
                    Cells(intRow + 1, 7) = _
                    WorksheetFunction.TTest(Range(Cells(intRegimeStart + 1, 2), Cells(intRow, 2)), _
                    Range(Cells(intRow + 1, 2), Cells(intRow + Cells(intRow + 1, 6), 2)), 2, 3)
                Else
                    Cells(intRow + 1, 7) = EqP(intRegimeStart, intRow)
                End If
            Else
                Cells(intRow + 1, 7) = ""
            End If
            intRegimeStart = intRow
        End If
    Next intRow
    'print the legend
    Cells(41, 10) = "RSI: Regime Shift Index"
    Cells(42, 10) = "Mean: Equal-weighed arithmetic means of the regimes"
    Cells(43, 10) = "Weighed: Weighed means of the regimes using the Huber's weight function " + _
    "with the parameter = " + frmMain.txtHuberParam.Text
    Cells(44, 10) = "Length: Length of the regimes"
                 If frmMain.optRNnone.Value = True Or _
                (frmMain.chkPrewhitening.Value = True And frmMain.chkFilteredData.Value = True) Then
        Cells(45, 10) = "Conf: Confidence level of the difference between the mean values of the " + _
        "neighboring regimes based on the Student's two-tailed t-test with unequal variance" + _
        " (TTEST procedure in Excel)."
    Else
        Cells(45, 10) = "Conf: Confidence level of the difference between the mean values of the " + _
        "neighboring regimes based on the Student's two-tailed t-test with unequal variance" + _
        " and equivalent sample size."
    End If
    Cells(46, 10) = "Outliers: Weight of the deviations from the weighed mean greater than " + _
    frmMain.txtHuberParam.Text + " standard deviation(s)."
    
    'plot the results
    Call RegimeChart(2)
    
   'print residuals
    wbkOut.Worksheets("Residuals").Activate
    
    'print labels in the first row
    If intTS = 2 Then
        ActiveSheet.UsedRange.Select
        Selection.Clear
    End If
    Range("a1").Select
    
    'print label in the first row
    If intTS = 2 Then
        Cells(1, 1) = Application.Workbooks(strActiveWbkName).Worksheets(strDataShtName).Cells(1, 1)
    End If
    Cells(1, intTS) = Application.Workbooks(strActiveWbkName).Worksheets(strDataShtName).Cells(1, intTS)
    
    'print time series
    If intTS = 2 Then
        For intRow = 1 To intMaxRows
            Cells(intRow + 1, 1) = rngMatrix(intRow + intDT - intLagF, 1)
        Next intRow
        Cells(intRow + 2, 1) = "Number of regimes"

    End If
    For intRow = 1 To N
        Cells(intRow + intStartRow, intTS) = sngResiduals(intRow)
    Next intRow
    Cells(intMaxRows + 3, intTS) = intNumberOfRegimes
    
    'plot the first time series of the residuals
    Call RegimeChart(3)
    
Summary_Sheet:
    wbkOut.Worksheets("Summary").Activate
    
    'print labels in the first row
    If intTS = 2 Then
        ActiveSheet.UsedRange.Select
        Selection.Clear
    End If
    Range("a1").Select
    
    'print label in the first row
    If intTS = 2 Then
        Cells(1, 1) = Application.Workbooks(strActiveWbkName).Worksheets(strDataShtName).Cells(1, 1)
    End If
    Cells(1, intTS) = Application.Workbooks(strActiveWbkName).Worksheets(strDataShtName).Cells(1, intTS)
    
    'print time series
    If intTS = 2 Then
        For intRow = 1 To intMaxRows
            Cells(intRow + 1, 1) = rngMatrix(intRow + intDT - intLagF, 1)
        Next intRow
    End If
    For intRow = 1 To N
        Cells(intRow + intStartRow, intTS) = Abs(rsi(intRow))
    Next intRow
    
Next intTS

'calculate summary rsi in the summary sheet
Cells(1, intMaxCols + 1) = "Sum"
Cells(1, intMaxCols + 2) = "Count"
For intRow = 2 To intMaxRows + 1 - intDT + intLagF
    Set rngW = Range(Cells(intRow, 2), Cells(intRow, intMaxCols))
    If WorksheetFunction.Count(rngW) > 0 Then
       Cells(intRow, intMaxCols + 1) = WorksheetFunction.sum(rngW)
    Else
       Cells(intRow, intMaxCols + 1) = ""
    End If
    'count shifts in a row
    intCountShifts = 0
    For intCol = 1 To intMaxCols - 1
        If rngW(intCol) > 0 Then intCountShifts = intCountShifts + 1
    Next intCol
    Cells(intRow, intMaxCols + 2) = intCountShifts
Next intRow
'plot Summary
Call SummaryChart(2)
Range("a1").Select

End Sub

Sub ShiftMean()
'Calculates rsi for the time series dblX

Dim n1 As Integer 'variable number of observations in regime 1
Dim dblR1() As Double
Dim dblDF As Double 'Degrees of freedom
Dim i As Integer 'counter
Dim dblVarVre() As Double 'holds subsample var

prob = frmMain.txtSignifLevel.Value
cutoff = frmMain.txtCutOffLen.Value

ReDim rsi(1 To N)
ReDim dblR1(1 To cutoff)
ReDim dblVarVre(1 To N - cutoff + 1)

'calculate average variance for a cutoff period
For intYear = cutoff To N
    For i = 1 To cutoff
        dblR1(i) = dblX(intYear - cutoff + i)
    Next i
    dblVarVre(intYear - cutoff + 1) = WorksheetFunction.Var(dblR1)
Next intYear
'dblAveVar = WorksheetFunction.Median(dblVarVre)
dblAveVar = WorksheetFunction.Average(dblVarVre)

'initial mean and diff
If frmMain.optRNnone = True Then
    dblDF = cutoff + cutoff - 2
Else
    If frmMain.chkPrewhitening.Value = True Then
        dblDF = cutoff + cutoff - 2
    Else
        dblDF = EqN(cutoff + cutoff - 2)
    End If
End If
dblDiff = WorksheetFunction.TInv(prob, dblDF) * Sqr(2 * dblAveVar / cutoff)
For i = 1 To cutoff
    dblR1(i) = dblX(i)
Next i
dblRegimeMean = WeighedAverage(dblR1)
ChangePoint = 1

'main cycle
For intYear = 2 To N
    If dblX(intYear) > dblRegimeMean + dblDiff Then
        rsi(intYear) = cusumUp
    Else
        If dblX(intYear) < dblRegimeMean - dblDiff Then
            rsi(intYear) = cusumDown
        Else
            rsi(intYear) = 0
        End If
    End If
    
    'check for the situation when the test is not over for the last
    'change point, but we are too close to the end of the time series
    If Abs(rsi(intYear)) > 0 And intYear > (N - cutoff + 1) Then Exit Sub
        
    If rsi(intYear) = 0 Then
    ' intYear is not a new changepoint
        If intYear >= ChangePoint + cutoff Then
            'recalculate regime mean and dblDiff
            'currently dblDiff remains constant for the entire process /series
            n1 = intYear - ChangePoint + 1
            ReDim dblR1(1 To n1)
            For i = 1 To n1
                dblR1(i) = dblX(ChangePoint + i - 1)
            Next i
            dblRegimeMean = WeighedAverage(dblR1)
            'old variant of dblDF:
            'dblDF = intYear - ChangePoint + cutoff - 1
            
            'now decision: effective N or not?
            If frmMain.optRNnone = True Then
                dblDF = cutoff + cutoff - 2
            Else
                If frmMain.chkPrewhitening.Value = True Then
                    dblDF = cutoff + cutoff - 2
                Else
                    dblDF = EqN(cutoff + cutoff - 2)
                End If
            End If
            dblDiff = WorksheetFunction.TInv(prob, dblDF) * Sqr(2 * dblAveVar / cutoff)
        End If
    Else
        'regime shift is detected
        ChangePoint = intYear
        'recalculate regime mean and dblDiff
        'currently dblDiff remains constant for the entire process /series
        ReDim dblR1(1 To cutoff)
        For i = 1 To cutoff
            dblR1(i) = dblX(ChangePoint + i - 1)
        Next i
        dblRegimeMean = WeighedAverage(dblR1)
        'decision: effective N?
        If frmMain.optRNnone = True Then
            dblDF = cutoff + cutoff - 2
        Else
            If frmMain.chkPrewhitening.Value = True Then
                dblDF = cutoff + cutoff - 2
            Else
                dblDF = EqN(cutoff + cutoff - 2)
            End If
        End If
        dblDiff = WorksheetFunction.TInv(prob, dblDF) * Sqr(2 * dblAveVar / cutoff)
    End If
Next intYear

End Sub
Function cusumUp() As Double
'This function calculates cusum for upward shifts
Dim i As Integer
Dim dblXdev As Double 'normalized deviation of x value from the expected value of the new
                        'regime, i.e., from dblRegimeMean + dblDiff
Dim dblXweight As Double 'Huber's weight of value x.
Dim dblSumofWeights As Double
                        
cusumUp = 0
For i = intYear To intYear + cutoff - 1
    If i > N Then
        If dblSumofWeights > 0 Then
            cusumUp = cusumUp / dblSumofWeights
            Exit Function
        End If
    End If
    
    dblXdev = (dblX(i) - dblRegimeMean - dblDiff) / Sqr(dblAveVar)
    
    'determine the weight of the normalized deviation
    If dblXdev = 0 Then
        dblXweight = 1
    Else
        dblXweight = WorksheetFunction.Min(1, frmMain.txtHuberParam.Value / Abs(dblXdev))
    End If
    
    'sum weights and weighed values
    dblSumofWeights = dblSumofWeights + dblXweight
    cusumUp = cusumUp + dblXdev * dblXweight
    
    'check if cusum turns zero
        If cusumUp < 0 Then
        cusumUp = 0
        Exit Function
    End If
Next i
cusumUp = cusumUp / dblSumofWeights

End Function
Function cusumDown() As Double
'This function calculates cusum fo downward shifts
Dim i As Integer
Dim dblXdev As Double 'normalized deviation of x value from the expected value of the new
                        'regime, i.e., from dblRegimeMean + dblDiff
Dim dblXweight As Double 'Huber's weight of value x.
Dim dblSumofWeights As Double

cusumDown = 0
For i = intYear To intYear + cutoff - 1
    If i > N Then
        If dblSumofWeights > 0 Then
            cusumDown = cusumDown / dblSumofWeights
            Exit Function
        End If
    End If
    
        dblXdev = (dblX(i) - dblRegimeMean + dblDiff) / Sqr(dblAveVar)
    
    'determine the weight of the normalized deviation
    If dblXdev = 0 Then
        dblXweight = 1
    Else
        dblXweight = WorksheetFunction.Min(1, frmMain.txtHuberParam.Value / Abs(dblXdev))
    End If
    
    'sum weights and weighed values
    dblSumofWeights = dblSumofWeights + dblXweight
    cusumDown = cusumDown + dblXdev * dblXweight
    
    'check if cusum turns zero
    If cusumDown > 0 Then
        cusumDown = 0
        Exit Function
    End If
Next i
cusumDown = cusumDown / dblSumofWeights

End Function
Function WeighedAverage(dblW() As Double) As Double
'Calculation of the mean estimate for a given range using Huber's weights
'The Huber's function is defined as min( 1, param/(|res|/scale))
'param: a given parameter (txtHuberParam); affects the range where weights = 1
'res: deviation from the mean estimate (I used the median as the first approximation)
'scale: estimate of variation. Some use 1.483*(median absolute deviation, or MAD,
' of the deviations of the data from their median). I used std, i.e. dblAveVar
'See http://www.stat.berkeley.edu/~stark/Preprints/Oersted/writeup.htm
'The procedure consists of two iterations. At the first iteration the estimate of the
'regime average is the simple unweighed arithmetic mean. In the second iteration
'it is weighed average from the first iteration.

Dim dblEstAve As Double 'estimated average of the range
Dim dblXdev As Double 'deviation from the median
Dim dblXweight As Double 'weights
Dim dblSumofWeights As Double
Dim dblSumAve As Double
Dim i As Integer
Dim k As Integer

dblEstAve = WorksheetFunction.Average(dblW)

For k = 1 To 2
  dblSumofWeights = 0
  dblSumAve = 0
  For i = 1 To UBound(dblW)
    dblXdev = (dblW(i) - dblEstAve) / Sqr(dblAveVar)
    'determine the weight of the normalized deviation
    If dblXdev = 0 Then
        dblXweight = 1
    Else
        dblXweight = WorksheetFunction.Min(1, frmMain.txtHuberParam.Value / Abs(dblXdev))
    End If
    
    'sum weights and weighed values
    dblSumofWeights = dblSumofWeights + dblXweight
    dblSumAve = dblSumAve + dblXdev * dblXweight
  Next i
  dblSumAve = dblSumAve / dblSumofWeights
  dblSumAve = dblSumAve * Sqr(dblAveVar) + dblEstAve
  dblEstAve = dblSumAve
Next k
WeighedAverage = dblEstAve

End Function
Public Function AlphaEst(x() As Double) As Double
'Calculate the bias-corrected estimate of alpha
'in x(t)= alpha * x(t-1) + white noise as the following:
' 1. Calculate est1, the OLS estimate of alpha for all
'   sequential subsamples of size nsub;
' 2. Let est1 be the median value of all est1 for subsamples
' 3. Make bias correction for est1

Dim ss() As Double 'to remember all est1 estimates for subsamples
Dim i As Integer 'counter
Dim j As Integer 'counter
Dim k As Integer 'counter
Dim y() As Double
Dim z() As Double
Dim est1 As Double ' OLS estimate of alpha


Nsub = frmMain.txtSubsampleSize.Value

'test the subsample size (Nsub) limits
If Nsub < 5 And frmMain.optRNKendall.Value = True Then
    Nsub = 5
    MsgBox "The subsample size is too small and changed to 5.", _
    vbInformation, "Incorrect subsample size"
End If

If Nsub < 3 And (frmMain.optIPN4.Value = True Or frmMain.optOLS.Value = True) Then
    Nsub = 3
    MsgBox "The subsample size is too small and changed to 3.", _
    vbInformation, "Incorrect subsample size"
End If

If Nsub > N Then Nsub = N

'Subsampling
ReDim ss(1 To N - Nsub + 1)
ReDim y(1 To Nsub)
For i = 1 To N - Nsub + 1
    For k = 1 To Nsub
        y(k) = x(i + k - 1)
    Next k
    ss(i) = OLSAR1(y)
Next i
est1 = WorksheetFunction.Median(ss)

' bias correction
If frmMain.optRNKendall = True Then
    AlphaEst = MPK(est1)
ElseIf frmMain.optIPN4 = True Then
    AlphaEst = IPN4(est1)
Else
    AlphaEst = est1
End If

'limits
If AlphaEst < 0 Then AlphaEst = 0
If AlphaEst > 0.99 Then AlphaEst = 0.99

End Function
Function OLSAR1(x() As Double) As Double
'OLS estimate of AR1 coefficient
'from Orcutt and Winokur, 1969, Econometrica, 37:1,1-14

Dim ave1 As Double 'first average from 2 to n
Dim ave2 As Double 'second
Dim sumNom As Double
Dim sumDenom As Double
Dim i As Integer
Dim Nobs As Integer 'number of observations in x()

Nobs = UBound(x)

For i = 2 To Nobs
    ave1 = ave1 + x(i)
    ave2 = ave2 + x(i - 1)
Next i
ave1 = ave1 / (Nobs - 1)
ave2 = ave2 / (Nobs - 1)
For i = 2 To Nobs
    sumNom = sumNom + (x(i) - ave1) * (x(i - 1) - ave2)
    sumDenom = sumDenom + (x(i - 1) - ave2) * (x(i - 1) - ave2)
Next i
If sumDenom > 0 Then
    OLSAR1 = sumNom / sumDenom
Else
    OLSAR1 = 0
End If

End Function
Function IPN4(est As Double) As Double
'Calculates bias-corrected AR(1) coefficient [est] using
'interative 4-step procedure with corrections that are
'reverse proportional to the subsample size - nsub
Dim i As Integer

IPN4 = est + 1 / Nsub
For i = 1 To 3
   IPN4 = IPN4 + Abs(IPN4) / Nsub
Next i

End Function

Function MPK(est As Double) As Double
'Calculates bias-corrected AR(1) coefficient using the formula of
'Marriott-Pope and Kendall(see Orcutt and Winokur, 1969, Econometrica, 37:1,1-14)
'est is the OLS estimate of AR(1)
'Nsub - sample size declared globaly

If Nsub > 4 Then
    MPK = ((Nsub - 1) * est + 1) / (Nsub - 4)
Else
    'should not be here
    MPK = est
End If

End Function
Function EqP(intRegimeStart As Integer, intRow As Integer) As Double
'calculates t-test using the equivalent sample size
'as in von Storch and Zwiers (1999, p. 115)
'Note: The formula for the equiv. sample size published in Zwiers and Storch (1995)
'in J. Climate is somewhat different

Dim rng1 As Range
Dim rng2 As Range
Dim ave1 As Double
Dim ave2 As Double
Dim var1 As Double
Dim var2 As Double
Dim eN1 As Integer 'effective sample size
Dim eN2 As Integer
Dim Ns As Integer 'sample size
Dim T_stat As Double 't-statistic
Dim i As Integer
Dim str As String

'define sample 1 and sample 2
Set rng1 = Range(Cells(intRegimeStart + 1, 2), Cells(intRow, 2))
Set rng2 = Range(Cells(intRow + 1, 2), Cells(intRow + Cells(intRow + 1, 6), 2))

'calculate the means and variances
ave1 = WorksheetFunction.Average(rng1)
ave2 = WorksheetFunction.Average(rng2)
var1 = WorksheetFunction.Var(rng1)
var2 = WorksheetFunction.Var(rng2)

'calculate effective sample sizes
'range1
Ns = rng1.Rows.Count
If Ns < 2 Then
    str = MsgBox("There is no data to calculate the effective sample size", vbCritical, _
    "Range not specified")
    EqP = -1
    Exit Function
End If
eN1 = EqN(Ns)

'range 2
Ns = rng2.Rows.Count
If Ns < 2 Then
    str = MsgBox("There is no data to calculate the effective sample size", vbCritical, _
    "Range not specified")
    EqP = -1
    Exit Function
End If
eN2 = EqN(Ns)

't-statistics
T_stat = Sqr(var1 / eN1 + var2 / eN2)
T_stat = Abs(ave1 - ave2) / T_stat

'probability for the t-distribution
EqP = WorksheetFunction.TDist(T_stat, eN1 + eN2 - 2, 2)

End Function
Function EqN(Ns As Integer) As Integer
'calculates the equivalent sample size
'as in von Storch and Zwiers (1999, p. 115)
'Note: The formula for the equiv. sample size published by
'Zwiers and Storch (1995)in J. Climate is somewhat different
'Ns original sample size
'alpha - AR1 coefficient defined for the entire module

Dim i As Integer
Dim sum As Double
For i = 1 To Ns - 1
    sum = sum + (1 - i / Ns) * alpha ^ i
Next i
EqN = Ns / (1 + sum)
'just in case
If EqN <= 2 Then EqN = 2
If EqN > Ns Then EqN = Ns

End Function

'''

if __name__ == "__main__":
  data=np.genfromtxt(fn,delimiter=",",names=True,filling_values =np.NaN) #fill with NaNs so must use their bounds
  x=data["B3"][:]       ####<<<<<<<<<<<<<<<<<<<<<<<<
  print AlphaEst(x, 20, option="optIPN4", returnmsgs=False)