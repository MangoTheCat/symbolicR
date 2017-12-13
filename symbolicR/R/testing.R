# git revision: $Rev: $
# Date of last change: $LastChangedDate: 26/06/2012 $
# Last changed by: $LastChangedBy: ccampbell $
# 
# Original author: ccampbell
# Copyright Mango Solutions, Chippenham, UK
###############################################################################

#' Executes the RUnit unit tests suite
#' 
#' contains unit tests distributed with the package
#' 
#' @param testPath A single string holding the path to the internal test suite scripts
#' @param protocolFile Single string with the name of the HTML report file produced by RUnit for the internal unit tests
#' @return standard internal test suite (object of class RUnitTestData returned by \code{runTestSuite} from the \code{RUnit} package)
#' @author Mango Solutions
#' @keywords debugging programming 
#' @examples 
#' \dontrun{ 
#'      runDDMoReMDLTests
#' }

runUnit.symbolicR <- function(testPath = system.file(package = "symbolicR", "testing"), protocolFile = "symbolicRtestResults.html")
{
    require(RUnit)
    
    # testSuite : object of class RUnitTestSTuite.  
    # This will describe the internal unit test suite to be executed
    testSuite <- defineTestSuite( name = "symbolicR Unit Test Suite", dirs = testPath, 
                                  testFileRegexp = "^runit.+\\.[rR]$",  testFuncRegexp = "^test.+" )
    
    # suiteResults : object of class RUnitTestData
    suiteResults <-  runTestSuite( testSuite )
    printHTMLProtocol( suiteResults, protocolFile )
    
    return(suiteResults)
}
