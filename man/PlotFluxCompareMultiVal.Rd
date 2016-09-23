\name{PlotFluxCompareMultiVal}
\alias{PlotFluxCompareMultiVal}
\title{Plot multiple time series comparing modeled and observed fluxes}
\usage{
PlotFluxCompareMultiVal(strDf.obs, strCol.obs = "q_cms", strDf.mod1,
  strCol.mod1 = "q_cms", strDf.mod2 = NULL, strCol.mod2 = "q_cms",
  strDf.mod3 = NULL, strCol.mod3 = "q_cms", strDf.mod4 = NULL, 
  strCol.mod4 = "q_cms", strDf.mod5 = NULL, strCol.mod5 = "q_cms",
  strDf.mod6 = NULL, strCol.mod6 = "q_cms",
  stdate = NULL, enddate = NULL, logy = FALSE, labelObs = "Observed",
  labelMod1 = "Model 1", labelMod2 = "Model 2",labelMod3 = "Model 3", 
  labelMod4 = "Model 4", labelMod5 = "Model 5", labelMod6 = "Model 6",
  title = "Observed and Modelled Fluxes")
}
\arguments{
\item{strDf.obs}{The OBSERVED flux time series dataframe (e.g., output from
\code{\link{ReadUsgsGage}}). The dataframe must contain a column of flux
values and a POSIXct column.}

\item{strCol.obs}{The name of the column containing the flux values for the
OBSERVED dataframe (DEFAULT="q_cms").}

\item{strDf.mod1}{The FIRST MODEL flux time series dataframe (e.g., output
from \code{\link{ReadFrxstPts}}). The dataframe must contain a column of
flux values and a POSIXct column.}

\item{strCol.mod1}{The name of the column containing the FIRST MODEL flux
values (DEFAULT="q_cms").}

\item{strDf.mod2}{The SECOND MODEL flux time series dataframe (e.g., output
from \code{\link{ReadFrxstPts}}). The dataframe must contain a column of
flux values and a POSIXct column.}

\item{strCol.mod2}{The name of the column containing the SECOND MODEL flux
values (DEFAULT="q_cms").}

\item{strDf.mod3}{The THIRD MODEL flux time series dataframe (e.g., output
from \code{\link{ReadFrxstPts}}). The dataframe must contain a column of
flux values and a POSIXct column.}

\item{strCol.mod3}{The name of the column containing the THIRD MODEL flux
values (DEFAULT="q_cms").}

\item{strDf.mod4}{The FOURTH MODEL flux time series dataframe (e.g., output
from \code{\link{ReadFrxstPts}}). The dataframe must contain a column of
flux values and a POSIXct column.}

\item{strCol.mod4}{The name of the column containing the FOURTH MODEL flux
values (DEFAULT="q_cms").}

\item{strDf.mod5}{The FIFTH MODEL flux time series dataframe (e.g., output
from \code{\link{ReadFrxstPts}}). The dataframe must contain a column of
flux values and a POSIXct column.}

\item{strCol.mod5}{The name of the column containing the FIFTH MODEL flux
values (DEFAULT="q_cms").}

\item{strDf.mod6}{The SIXTH MODEL flux time series dataframe (e.g., output
from \code{\link{ReadFrxstPts}}). The dataframe must contain a column of
flux values and a POSIXct column.}

\item{strCol.mod6}{The name of the column containing the SIXTH MODEL flux
values (DEFAULT="q_cms").}

\item{stdate}{Start date for plot/statistics (DEFAULT=NULL, all records will
be used). Date MUST be specified in POSIXct format with appropriate
timezone (e.g., as.POSIXct("2013-05-01 00:00:00", format="\%Y-\%m-\%d
\%H:\%M:\%S", tz="UTC"))}

\item{enddate}{End date for plot/statistics (DEFAULT=NULL, all records will
be used). Date MUST be specified in POSIXct format with appropriate
timezone (e.g., as.POSIXct("2013-05-01 00:00:00", format="\%Y-\%m-\%d
\%H:\%M:\%S", tz="UTC"))}

\item{logy}{(TRUE or FALSE) Optional flag to set the y-axis to log-scale
(DEFAULT=FALSE).}

\item{labelObs}{Optional label for the observed streamflow
(DEFAULT="Observed")}

\item{labelMod1}{Optional label for the FIRST MODEL (DEFAULT="Model 1")}

\item{labelMod2}{Optional label for the SECOND MODEL (DEFAULT="Model 2")}

\item{labelMod3}{Optional label for the THIRD MODEL (DEFAULT="Model 3")}

\item{labelMod4}{Optional label for the FOURTH MODEL (DEFAULT="Model 4")}

\item{labelMod5}{Optional label for the FIFTH MODEL (DEFAULT="Model 5")}

\item{labelMod6}{Optional label for the SIXTH MODEL (DEFAULT="Model 6")}

\item{title}{Optional for the plot (DEFAULT="Observed and Modelled Fluxes")}
}
\value{
A plot of the hydrographs.
}
\description{
\code{PlotFluxCompare} plots a time series of an observed flux (e.g.,
streamflow, ET) and up to 6 modelled fluxes.
}
\details{
\code{PlotFluxCompare} reads modelled and observed dataframes (e.g., as
generated from \code{\link{ReadFrxstPts}} and \code{\link{ReadUsgsGage}}) and
plot the time series and summary statistics. The tool will subset data to
matching time periods (e.g., if the observed data is at 5-min increments and
modelled data is at 1-hr increments, the tool will subset the observed data
to select only observations on the matching hour break).
}
\examples{
## Take a time series of observed 5-minute streamflow values for Fourmile
## Creek (obsStr5min.fc) and two model runs (mod1Str1h.fc, mod2Str1h.fc),
## all with streamflow columns named "q_cms", and plot the the hydrographs
## for all three over the May-June snowmelt period.

\dontrun{
PlotFluxCompareMultiVal(obsStr5min.fc, "q_cms", modStrh.chrt.fc, "q_cms",
                strDf.mod2=modStrh.allrt.fc, strCol.mod2="q_cms",
                labelObs="Observed Fourmile Creek at Orodell",
                labelMod1="Channel Routing Only", labelMod2="All Routing",
                title="Streamflow: Fourmile Creek",
                stdate=as.POSIXct("2013-05-01 00:00:00",
                                  format="\%Y-\%m-\%d \%H:\%M:\%S", tz="UTC"),
                enddate=as.POSIXct("2013-06-30 00:00:00",
                                   format="\%Y-\%m-\%d \%H:\%M:\%S", tz="UTC"))
}
}
