# powershell script
#
# Run the HGS tool chain (grok, phgs, hsplot), plus postprocessing* sripts or
# preplot
#

$console = "console.txt"
$prefixfile = "batch.pfx"
$scriptForMerge = "" # skip merge with 0-string
#$scriptForMerge = "tp-mergeHSPlotOutputs.py"
$preprocessFilePattern  = "preprocess*"
$postprocessFilePattern = "postprocess*"

$EXE_GROK = "grok.exe"
$EXE_HGS = "phgs.exe"
$EXE_HSPLOT = "hsplot.exe"

########################################
# Helper functions
function Run-Processing($pattern) {
   $code=0
   # sort scripts
   #    e.g. postprocess.0.ps1 runs before postprocess.1.py
   $ssc = ( Get-ChildItem $pattern | sort )
   foreach ( $sc in $ssc ) {
      # check executable status
      if( get-command $sc -ErrorAction SilentlyContinue ) {
         $short_sc="$((get-item $sc).Name)"
         echo "--- Running $short_sc ---"

         if ( $sc -match '.py$' ) {
            # prepend 'python' so the correct interpreter is launched
            # (useful for virtual environments)
            & python $sc >> $console
         }
         else {
            & $sc >> $console
         }

         if( $LastExitCode -eq 0 ) {
            echo "--- Normal exit: $short_sc ---"
         }
         else {
            echo "--- Failed: $short_sc ---"
            if ( $ssc.length -gt 1 ) {
               echo "--- stopping $($pattern) script sequence ==="
            }
            $code=1
            break
         }
      }
      else {
         echo "--- not running $sc - not set as executable ---"
      }
   }
   $LastExitCode = $code
}

########################################
# Warn.
echo ""
echo "Disable Properties->QuickEdit to avoid console (script execution) hangs."
echo ""
# TODO
# http://stackoverflow.com/questions/30872345/script-commands-to-disable-quick-edit-mode

########################################
# get the problem prefix
$prob=""
if ( Test-Path $prefixfile ) {
   # assume the prefix is stored in this file
   $prob= Get-Content $prefixfile
}
else {
   # prompt and write prefix for next time
   $prob= Read-Host -Prompt 'Input problem prefix'
   echo "Writing prefix file $prefixfile" ""
   echo $prob | Out-File -FilePath $prefixfile -Encoding ascii
}

########################################
# print the start time
echo "`n--- Starting hgs-runall $(date) ---`n"
echo "`n--- Starting hgs-runall $(date) ---`n" > $console

########################################
# run preprocessing
if( Test-Path $preprocessFilePattern ) {
   echo "`n--- Starting preprocessing $(date) ---`n"
   Run-Processing $preprocessFilePattern

   if( $LastExitCode -eq 0 ) {
   echo "--- preprocessing done $(date) ---"
   }
   else {
      echo "--- prepocessing failed $(date) ---"
      exit 1
   }
}



########################################
# run grok
echo "`n--- Starting grok $(date) ---`n"

if ( ( & $EXE_GROK | Select-Object -Last 3 | Select-Object -First 1 )  -notmatch "normal exit" ) {
   echo "--- grok failed, not running phgs ---"
   echo ( Get-Content ($prob+'o.eco') | Select-Object -Last 10 )
   date
   exit 1
}
else {
   echo "--- grok: normal exit ---" ""
   [Console]::Out.Flush()
}


########################################
# run phgs
echo "`n--- Starting phgs $(date) ---`n"

# Works, but writes a big 'console' file
#echo "phgs console output is being sent to $console" ""
#phgs > $console
#$phgsstatus = Get-Content -Tail 17 $console

# TODO: print output to console AND get the last 17 lines
# TODO: write console output to a buffer
#$phgsstatus = Get-Content -Tail 17 phgs
#phgs | Tee-Object | Select-Object -Last 17

# test this... does it work?
echo "phgs console output suppressed. See $($prob)o.lst for progress"
echo ""

$phgsstatus = ( ( & $EXE_HGS ) | select-object -last 17 )

echo $phgsstatus[0..17]
if ( $phgsstatus | Select-String -Pattern "NORMAL EXIT" ) {
   echo "--- phgs: normal exit ---" ""
   [Console]::Out.Flush()
}
else {
   echo "--- phgs failed; not running hsplot ---"
   date
   exit 1
}


########################################
# run hsplot
#
echo "`n--- Starting hsplot $(date) ---`n"

if ( ( & $EXE_HSPLOT | Select-Object -Last 1 )  -notmatch "Normal Exit" ) {
   echo "--- hsplot failed, not running preplot(tecplot) ---"
   date
   exit 1
}
else {
   echo "--- hsplot: normal exit ---" ""
   [Console]::Out.Flush()
}

date

#
# get ready to generate tecplot-ready data files
#

# A list of files that will be sent to preplot in this script
# If other postprocessing scripts exist, this will be ignored
$preplotInList = @($prob + 'o.pm.dat')

if( test-path *o.dual.dat ) {
      $preplotInList += $prob+"o.dual.dat"
}

########################################
# merge the two output files
#if( ($scriptForMerge) -and (test-path *o.frac.dat) ) {
#
#   if( get-command $scriptForMerge -ErrorAction SilentlyContinue ) {
#
#      echo "`n--- Starting $scriptForMerge $(date) ---`n"
#         
#      $tmpf = [System.IO.Path]::GetTempFileName()
#      & $scriptForMerge --of $tmpf
#      $preplotInList = @( $prob+'o.dat' )
#
#      if( $LastExitCode -eq 0 ) {
#         echo "" "--- normal exit, $scriptForMerge ---" ""
#         if( test-path $preplotInList[0] ) { rm $preplotInList[0] }
#         mv $tmpf $preplotInList[0]
#      }
#      else {
#         echo "" "--- failed, $scriptForMerge ---" ""
#         rm $tmpf
#         if( test-path $preplotInList[0] ) { echo "Warning: keeping old $preplotInList[0]" }
#      }
#   }
#   else {
#      $preplotInList += $prob+"o.frac.dat"
#   }
#}
#
#date

########################################
# run scripts or run preplot
#
#
# look for and execute postprocess scripts
# else, run preplot on the spacial output only
#

if( Test-Path $postprocessFilePattern ) {
   echo "`n--- Starting postprocessing $(date) ---`n"
   Run-Processing $postprocessFilePattern
   echo "`n--- postprocessing done $(date) ---`n"
}
#else {
#   echo "`n--- Starting preplot $(date) ---`n"
#
#   foreach ( $plotin in $preplotInList ) {
#      $plotout = $plotin.substring(0,$plotin.length-3)+"plt"
#      # run preplot to prepare tecplot binary files
#      echo "--- running preplot $plotin $plotout >> $console ---"
#      & preplot $plotin $plotout >> $console
#
#      if( $LastExitCode -eq 0 ) {
#         echo "" "--- preplot: normal exit. Tecplot binary data is in $plotout. ---" ""
#      }
#      else {
#         echo "--- preplot failed. See $console. ---"
#      }
#   }
#}


########################################
# Clean up annoying files
# move all <prefix>o.*.NNNN files
# remove scratch_ files
#echo "`n--- Starting output file re-organization $(date) ---`n"
#mkdir -Force "$($prob)o.NNNN" | Out-Null
#mv -Force "$($prob)o.*.[0-9][0-9][0-9][0-9]" "$($prob)o.NNNN"
del scratch_*
