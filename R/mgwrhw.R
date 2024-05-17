#' mgwrhw
#'
#' displays the GWR and mixed GWR models automatically along with the tests and significance maps that are formed.
#'
#' @param dpk dataframe all variables that come from the shp data format and have geometric attributes that are usually imported with the st_read function from library(sf)
#' @param pers.reg The form of the regression equation that will be used as a GWR model is in the general form y~x1+x2+x3
#' @param coor_lat the name of the variable that is in the dpk dataframe that contains latitude coordinates and is written with quotation marks such as "Latitude" which indicates a column named Latitude
#' @param coor_long the name of the variable that is in the dpk dataframe that contains latitude coordinates and is written with quotation marks such as "Longitude" which indicates a column named Longitude
#' @param vardep the name of a variable that is in a dpk dataframe that contains one dependent variable and is written with quotation marks such as "y" which indicates a column named y
#' @param GWRonly user option to choose to display GWR results only or to form an MGWR model. Option 1 displays GWR output only while option 0 displays GWR and MGWR output.
#' @param kp user option to select kernel functions. Option 1 for Fixed Bisquare, option 2 for Fixed Gaussian, option 3 for Adaptive Bisquare, and option 4 for Adaptive Bisquare
#' @param alp alpha value (type 1 error) used in spatial regression model
#'
#' @export
#'
#' @import spgwr ggplot2 dplyr tidyr sf
#' @importFrom stats anova lm pf pt qf qt
#' @importFrom psych tr
#' @importFrom utils capture.output
#'
#'
#' @returns no return value, called for side effects
#' @return This function returns a list with the following objects:
#'
#'
#' ## for Mixed GWR model (GWRonly = 0)
#'
#' the general equation form of the Mixed GWR model is
#'
#' \eqn{y_{i}} = \eqn{\beta_{0}}(\eqn{u_{i}},\eqn{v_{i}}) + \eqn{\sum}\eqn{\beta_{k}}(\eqn{u_{i}},\eqn{v_{i}})\eqn{x_{ik}} + \eqn{\sum}\eqn{\beta_{k}}\eqn{x_{ik}} + \eqn{\epsilon_{i}}
#'
#' \describe{
#'   \item{output}{A character vector containing the captured output of GWR model and Mixed GWR model.}
#'   \item{gwr}{The result of the GWR model include CV, bandwith, Quasi R square, etc.}
#'   \item{Variability.Test}{Results of the variability test for global and local variables.}
#'
#' \eqn{H_{0}} : \eqn{\beta_{k}}(\eqn{u_{1}},\eqn{v_{1}}) = \eqn{\beta_{k}}(\eqn{u_{2}},\eqn{v_{2}}) = ... = \eqn{\beta_{k}}(\eqn{u_{n}},\eqn{v_{n}})
#'
#' \eqn{H_{1}} : not all \eqn{\beta_{k}}(\eqn{u_{i}},\eqn{v_{i}}) (\eqn{i} = 1, 2, ..., n) are equal
#'
#'
#' \deqn{F_{Variability.Test_{k}} = \frac{V^{2}_{k}{/}\gamma_{1}}{\widehat{\sigma}}}
#'
#' Conclusion : Reject \eqn{H_{0}} if \eqn{F_{Variability.Test_{k}}} \eqn{\geq} \eqn{F_{\alpha}}(\eqn{\frac{\gamma_{1}^{2}}{\gamma_{2}},\frac{\delta_{1}^{2}}{\delta_{2}}}) or p-value < \eqn{\alpha}.
#'
#' If \eqn{H_{0}} is rejected, it means that the k-th variable has a local influence, while if \eqn{H_{0}} fails to be rejected, it means that the k-th variable has a global influence.
#'
#' Reference : Leung, Y., Mei, C.L., & Zhang, W.X., (2000). "Statistic Tests for Spatial Non-Stationarity Based on the Geographically Weighted Regression Model", Environment and Planning A, 32 pp. 9-32. <doi:10.1068/a3162>.
#'
#'
#'   \item{F1.F2.F3.mgwr.Test}{Results of the F1(GoF Mixed GWR), F2(Global Simultaneous), F3(Local Simultaneous) tests.}
#'
#'   F1(GoF Mixed GWR) :
#'
#'   \eqn{H_{0}} : \eqn{\beta_{k}}(\eqn{u_{i}},\eqn{v_{i}}) = \eqn{\beta_{k}}
#'
#'   \eqn{H_{1}} : at least there is one  \eqn{\beta_{k}}(\eqn{u_{i}},\eqn{v_{i}}) \eqn{\neq} \eqn{\beta_{k}}
#'
#'   \deqn{F(1) = \frac{y^{T}((I-H)-(I-S)^{T}(I-S))y {/} v_{1}} {y^{T}(I-S)^{T}(I-S)y {/} u_{1}}}
#'
#'   if \eqn{H_{0}} is rejected, it shows that the Mixed GWR model is different from the OLS model]
#'
#'   F2(Global Simultaneous) :
#'
#'   \eqn{H_{0}} : \eqn{\beta_{q+1}} = \eqn{\beta_{q+2}} = ... = \eqn{\beta_{p}} = 0
#'
#'   \eqn{H_{1}} : at least one of \eqn{\beta_{k}} \eqn{\neq} 0
#'
#'   \deqn{F(2) = \frac{y^{T}((I-S_{l})^{T}(I-S_{l})-(I-S)^{T}(I-S))y {/} r_{1}} {y^{T}(I-S)^{T}(I-S)y {/} u_{1}}}
#'
#'   If \eqn{H_{0}} is rejected, it indicates that there is at least one global variable that has a significant effect in the model
#'
#'   F3(Local Simultaneous)
#'
#'   \eqn{H_{0}} : \eqn{\beta_{1}}(\eqn{u_{i}},\eqn{v_{i}}) = \eqn{\beta_{2}}(\eqn{u_{i}},\eqn{v_{i}}) = ... = \eqn{\beta_{q}}(\eqn{u_{i}},\eqn{v_{i}}) = 0
#'
#'   \eqn{H_{1}} : at least one of \eqn{\beta_{k}}(\eqn{u_{i}},\eqn{v_{i}}) \eqn{\neq} 0
#'
#'   \deqn{F(2) = \frac{y^{T}((I-S_{g})^{T}(I-S_{g})-(I-S)^{T}(I-S))y {/} r_{1}} {y^{T}(I-S)^{T}(I-S)y {/} u_{1}}}
#'
#'   If \eqn{H_{0}} is rejected, it indicates that there is at least one local variable that has a significant effect in the model
#'
#'   Reference : Yasin, & Purhadi. (2012). "Mixed Geographically Weighted Regression Model (Case Study the Percentage of Poor Households in Mojokerto 2008)". European Journal of Scientific Research, 188-196.
#'   <https://www.researchgate.net/profile/Hasbi-Yasin-2/publication/289689583_Mixed_geographically_weighted_regression_model_case_study_The_percentage_of_poor_households_in_Mojokerto_2008/links/58e46aa40f7e9bbe9c94d641/Mixed-geographically-weighted-regression-model-case-study-The-percentage-of-poor-households-in-Mojokerto-2008.pdf>.
#'
#'
#'   \item{Global.Partial.Test}{Results of the global partial test.}
#'
#'   \eqn{H_{0}} : \eqn{\beta_{k}} = 0 (k-th global variables are not significant)
#'
#'   \eqn{H_{1}} : \eqn{\beta_{k}} \eqn{\neq} 0 (k-th global variables are significant)
#'
#'   \deqn{T_{g} = \frac{\widehat{\beta_{k}}}{\widehat{\sigma}\sqrt{g_{kk}}}}
#'
#'   If \eqn{H_{0}} is rejected, it indicates that the k-th global variable has a significant effect
#'
#'   Reference : Yasin, & Purhadi. (2012). "Mixed Geographically Weighted Regression Model (Case Study the Percentage of Poor Households in Mojokerto 2008)". European Journal of Scientific Research, 188-196.
#'   <https://www.researchgate.net/profile/Hasbi-Yasin-2/publication/289689583_Mixed_geographically_weighted_regression_model_case_study_The_percentage_of_poor_households_in_Mojokerto_2008/links/58e46aa40f7e9bbe9c94d641/Mixed-geographically-weighted-regression-model-case-study-The-percentage-of-poor-households-in-Mojokerto-2008.pdf>.
#'
#'
#'   \item{map.mgwr}{Visualization of Mixed GWR results in the form of a regional map with variables that are significant globally and locally.}
#'   \item{Global_variable}{A list of global variables used in the analysis.}
#'   \item{Local_variable}{A list of local variables used in the analysis.}
#'   \item{AICc}{The corrected Akaike Information Criterion.}
#'   \item{AIC}{The Akaike Information Criterion.}
#'   \item{R_square}{The coefficient of determination.}
#'   \item{adj_R_square}{The adjusted coefficient of determination.}
#'   \item{table.mgwr}{A data frame about output table of MGWR model (include estimator, standar error, t-statistics, p-value).}
#' }
#'
#' ## for GWR model (GWRonly = 1)
#'
#' the general equation form of the GWR model is
#'
#' \eqn{y_{i}} = \eqn{\beta_{0}}(\eqn{u_{i}},\eqn{v_{i}}) + \eqn{\sum}\eqn{\beta_{k}}(\eqn{u_{i}},\eqn{v_{i}})\eqn{x_{ik}} + \eqn{\epsilon_{i}}
#'
#' \describe{
#'   \item{output}{A character vector containing the captured output of GWR model.}
#'   \item{gwr}{A character vector containing the result of the GWR model include CV, bandwith, Quasi R square, etc.}
#'   \item{GoF.test}{A character vector containing the results of the Godness of Fit Test.}
#'   \item{anova_gwr}{Results of the anova table.}
#'   \item{map.gwr}{Visualization of the GWR results.}
#'   \item{table.gwr}{A data frame about output table of GWR model (include estimator, standar error, t-statistics, p-value).}
#' }
#'
#'
#' @examples
#' mod1 = mgwrhw(dpk=redsb, pers.reg = Y ~ X2 + X4 + X5 + X6,
#' coor_lat = "Latitude", coor_long = "Longitude",
#' vardep = "Y", GWRonly = 0, kp = 3, alp = 0.05)
#' mod1$gwr
#' mod1$Variability.Test
#' mod1$Global_variable
#' mod1$Local_variable
#' mod1$F1.F2.F3.mgwr.Test
#' mod1$Global.Partial.Test
#' mod1$map.mgwr

mgwrhw = function(dpk, pers.reg, coor_lat, coor_long, vardep, GWRonly, kp, alp){
  dtk = data.frame(dpk) #proses mengubah spatialdataframe ke dataframe biasa
  dt = dtk[,-length(dtk)] #mengahpus geometri agar tidak mengganggu saat memanggil satu variabel
  y=as.matrix(dt[vardep])
  lat=as.matrix(dt[coor_lat])
  lon=as.matrix(dt[coor_long])
  n=nrow(dt[vardep])
  I=diag(1,n,n)
  W=array(0,dim=c(n,n,n))
  d=matrix(0,n,n)
  reg.klasik=lm(pers.reg,data=dt)
  result.um = list()


  if (kp==1){
    # FIXED GAUSSIAN #############################################################################################################
    #bandwidth
    fixgauss= spgwr::gwr.sel(pers.reg,data=dt,coords=cbind(lon,lat),adapt=FALSE,gweight=gwr.Gauss)
    #estimasi parameter
    gwr1=spgwr::gwr(pers.reg,data=dt,coords=cbind(lon,lat),
             bandwidth=fixgauss,hatmatrix=TRUE,gweight=gwr.Gauss)

    bw=gwr1$bandwidth #menampilkan nilai bandwidth

    ##############
    ##### FG #####
    ##############

    for (i in 1:n)
    {
      for (j in 1:n)
      {
        d[i,j]=sqrt((lat[i,1]-lat[j,1])^2+(lon[i,1]-lon[j,1])^2)

        W[j,j,i]=exp(-1/2*((d[i,j]/bw)^2))
      }
    }

  } else {

    if(kp==2){
      # FIXED BISQUARE #############################################################################################################
      #bandwidth
      fixbisquare=spgwr::gwr.sel(pers.reg,data=dt,coords=cbind(lon,lat),adapt=FALSE,gweight=gwr.bisquare)
      #estimasi parameter
      gwr1=spgwr::gwr(pers.reg,data=dt,coords=cbind(lon,lat),
               bandwidth=fixbisquare,hatmatrix=TRUE,gweight=gwr.bisquare)

      bw=gwr1$bandwidth #menampilkan nilai bandwidth

      ##############
      ##### FB #####
      ##############

      for (i in 1:n)
      {
        for (j in 1:n)
        {
          d[i,j]=sqrt((lat[i,1]-lat[j,1])^2+(lon[i,1]-lon[j,1])^2)
          if(d[i,j] < bw){
            W[j,j,i]=(1-((d[i,j]/bw)^2))^2} else
            { W[j,j,i] = 0 }
        }
      }


    } else {

      if(kp==3){

        # ADAPTIVE GAUSSIAN #############################################################################################################
        #bandwidth
        adaptgauss=spgwr::gwr.sel(pers.reg,data=dt,coords=cbind(lon,lat),adapt=TRUE,gweight=gwr.Gauss)
        #estimasi parameter
        gwr1=spgwr::gwr(pers.reg,data=dt,coords=cbind(lon,lat),
                 adapt=adaptgauss,hatmatrix=TRUE,gweight=gwr.Gauss)

        bw=gwr1$bandwidth #menampilkan nilai bandwidth

        ##############
        ##### AG #####
        ##############

        for (i in 1:n)
        {
          for (j in 1:n)
          {
            d[i,j]=sqrt((lat[i,1]-lat[j,1])^2+(lon[i,1]-lon[j,1])^2)

            W[j,j,i]=exp(-1/2*((d[i,j]/bw[i])^2))
          }
        }


      } else {

        if(kp==4){

          # ADAPTIVE BISQUARE #############################################################################################################
          #bandwidth
          adaptbisquare=spgwr::gwr.sel(pers.reg,data=dt,coords=cbind(lon,lat),adapt=TRUE,gweight=gwr.bisquare)
          #estimasi parameter
          gwr1=spgwr::gwr(pers.reg,data=dt,coords=cbind(lon,lat),
                   adapt=adaptbisquare,hatmatrix=TRUE,gweight=gwr.bisquare)

          bw=gwr1$bandwidth #menampilkan nilai bandwidth

          ##############
          ##### AB #####
          ##############

          for (i in 1:n)
          {
            for (j in 1:n)
            {
              d[i,j]=sqrt((lat[i,1]-lat[j,1])^2+(lon[i,1]-lon[j,1])^2)
              if(d[i,j] < bw[i]){
                W[j,j,i]=(1-((d[i,j]/bw[i])^2))^2} else
                { W[j,j,i] = 0 }
            }
          }

        }

      }

    }

  }

  ## output model terpilih
  Fnya = spgwr::BFC02.gwr.test(gwr1) #H0 : model GWR = model OLS
  aan = anova(reg.klasik)

  #----------------------------------------------------------------------------------------
  qradj = 1 - ((Fnya$estimates[2]/sum(aan$`Sum Sq`)*(sum(aan$Df)/Fnya$parameter[2])))
  ### ----------------------------------------------------------------------- ###

  ##GRAFIK UNTUK GWR
  db_t = gwr1$results$edf  #sudah aman tidak perlu Snya lagi
  df_t_gwr = data.frame()

  for(i in 1:length(colnames(gwr1$lm$x))){
    for(j in 1:nrow(dt[vardep])){
      df_t_gwr[j,i+(i-1)]=data.frame(gwr1$SDF)[j,i+1]/data.frame(gwr1$SDF)[j,i+1+length(colnames(gwr1$lm$x))]

      if(df_t_gwr[j,i+(i-1)]>0){
        df_t_gwr[j,i+i]=2*pt(df_t_gwr[j,i+(i-1)],df=db_t,lower.tail=FALSE)} else
        {df_t_gwr[j,i+i]=2*pt((-1)*df_t_gwr[j,i+(i-1)],df=db_t,lower.tail=FALSE)}

    }
  }

  ## kasih nama df_t_gwr
  varig = c("t_","pv_")
  vdp = rep(varig,length(names(gwr1$lm$coefficients)))
  gdup = rep(names(gwr1$lm$coefficients),each=length(varig))
  nco=c()
  ggd = data.frame(vdp,gdup)%>% tidyr::unite(nco, c(1:2), sep = "")

  for(i in 1:length(colnames(gwr1$lm$x))){
    for(j in 1:nrow(dt[vardep])){
      if(df_t_gwr[j,i+i] < alp) {
        df_t_gwr[j,length(colnames(gwr1$lm$x))*2+i] = colnames(gwr1$lm$x)[i]
      } else {
        df_t_gwr[j,length(colnames(gwr1$lm$x))*2+i] = " "
      }

    }
  }

  vdmgw = rep("sig_",length(names(gwr1$lm$coefficients)))
  vdmugw = data.frame(vdmgw,names(gwr1$lm$coefficients))%>%tidyr::unite(nco, c(1:2),sep = "")
  df_t_gwr = `colnames<-`(df_t_gwr,c(ggd$nco,vdmugw$nco))

  GWRall = df_t_gwr %>%
    tidyr::unite(col = klasgwr, c((length(colnames(gwr1$lm$x))*2+2):(length(colnames(df_t_gwr)))), sep = "/")
  klasgwr = c()
  tablegwr = data.frame(dt, as.data.frame(gwr1$SDF), df_t_gwr, klasgwr=GWRall$klasgwr)

  VGWR = ggplot2::ggplot(dpk)+
    geom_sf(color = "white", aes(fill = GWRall$klasgwr)) +
    theme(legend.position = "left")


  print.mgwrhwc <- function(x, ...) {
    return(invisible(x))
  }

  if (GWRonly==1){
    z <- list(
      output = capture.output({
        cat("                                                 ","\n")
        cat("                                                 ","\n")
        cat("===================== Output GWR ================","\n")
        cat("                                                 ","\n")
        print(gwr1)
        cat("Adj. R squared = ",qradj,"\n")
        cat("                                                 ","\n")
        cat("Godness Of Fit Test:","\n")
        print(Fnya) #H0 : model GWR = model OLS
        cat("                                                 ","\n")
        cat("ANOVA:","\n")
        print(spgwr::anova.gwr(gwr1))
        cat("===================== ========== ================","\n")
        cat("                                                 ","\n")
      }),
      gwr=gwr1, GoF.test = Fnya,
      anova_gwr = spgwr::anova.gwr(gwr1),
      table.gwr = tablegwr,
      map.gwr = VGWR)

    class(z) <- "mgwrhwc"
    print.mgwrhwc(z)

}
   else {
    if(GWRonly==0){

      ## MGWR ###----------------------------------------------------


      #uji variabilitas untuk cari variabel global
      LF3 = spgwr::LMZ.F3GWR.test(gwr1) # yg gak signifikan itu ln_perkap_mak, kawinpr, imunisasi # cara menyembunyikan ini gimana???

      varlok = c()
      varglob = c()
      for(i in 2:length(reg.klasik$model)){
        if(LF3[i,4]<alp){
          varlok = c(varlok,names(LF3[,4])[i])

        } else {varglob = c(varglob,names(LF3[,4])[i])}
      }


      bw=gwr1$bandwidth
      intercept = rep(1,length(y))
      xl = as.matrix(cbind(intercept,dt[varlok])) #berhasil
      xg = as.matrix(cbind(dt[varglob]))
      x = as.matrix(cbind(xl,xg))

      ng=ncol(xg)
      nl=ncol(xl)
      n=length(y)

      ## Estimasi Parameter Global


      Sl=matrix(0,n,n) # definisi Sl

      for (i in 1:n)
      {
        Sl[i,] = xl[i,] %*% (solve(t(xl)%*%(W[,,i])%*%xl)) %*% t(xl) %*% (W[,,i])
      }

      beta.g=( (solve((t(xg)%*%t(I-Sl))%*%(I-Sl)%*%xg)) %*% t(xg) %*% t(I-Sl) %*% (I-Sl) ) %*% y


      ## Estimasi Parameter Lokal
      beta.l=matrix(0,n,nl) #MATRIKS BETA LOKAL
      for (i in 1:n)
      {
        beta.l[i,]=t( ( solve((t(xl)%*%(W[,,i]))%*%xl) %*% t(xl) %*% (W[,,i]) )%*%(y-(xg%*%beta.g)) )
      }

      S=Sl + ( (I-Sl) %*% xg %*% solve(((t(xg)%*%t(I-Sl))%*%(I-Sl))%*%xg) %*% t(xg) %*% t(I-Sl) %*%(I-Sl) )

      ## Prediksi
      y.hat=S%*%y
      #y.hat #PREDIKSI
      residual=(I-S)%*%y

      ## Pengujian

      H = x %*% solve(t(x)%*%x) %*% t(x)
      Sg = xg %*% solve(t(xg)%*%xg) %*% t(xg)


      v=c(0,0)
      u=c(0,0)
      r=c(0,0)
      t=c(0,0)
      for (i in 1:2)
      {
        v[i]=psych::tr(((I-H)-(t(I-S)%*%(I-S)))^i)
        u[i]=psych::tr((t(I-S)%*%(I-S))^i)
        r[i]=psych::tr( ((t(I-Sl)%*%(I-Sl)) - (t(I-S)%*%(I-S)))^i )
        t[i]=psych::tr( ((t(I-Sg)%*%(I-Sg)) - (t(I-S)%*%(I-S)))^i )
      }


      F1=as.vector( (((t(y)%*%((I-H)-(t(I-S)%*%(I-S))))%*%y)/v[1]) /
                      ((t(y)%*%t(I-S)%*%(I-S)%*%y)/u[1]) )
      df1.1=(v[1]^2/v[2])
      df2=(u[1]^2)/u[2]

      F2=as.vector( ( (t(y) %*% ((t(I-Sl)%*%(I-Sl)) - (t(I-S)%*%(I-S))) %*% y)/r[1] ) /
                      ((((t(y)%*%t(I-S))%*%(I-S))%*%y)/u[1]) )
      df1.2=(r[1]^2/r[2])

      F3=as.vector( ( (t(y) %*% ((t(I-Sg)%*%(I-Sg)) - (t(I-S)%*%(I-S))) %*% y)/t[1] ) /
                      ((((t(y)%*%t(I-S))%*%(I-S))%*%y)/u[1]) )
      df1.3=(t[1]^2/t[2])



      ### Pengujian GoF MGWR, Simultan Global dan Simultan Lokal

      F=as.vector(rbind(F1,F2,F3))
      df1=c(df1.1,df1.2,df1.3)
      f.tabel=as.vector(matrix(0,3,1))
      p.value=as.vector(matrix(0,3,1))
      for (i in 1:3)
      {
        f.tabel[i]=qf((1-alp), df1=df1[i], df2=df2)
        p.value[i]=1-(pf(F[i], df1=df1[i], df2=df2))
      }

      Uji.Serentak=cbind(F,f.tabel,df1,df2,p.value=round(p.value,11))


      ### Pengujian Parsial Global
      G = solve(t(xg)%*%t(I-Sl)%*%(I-Sl)%*%xg) %*% t(xg) %*% t(I-Sl) %*% (I-Sl)
      gkk = diag(G%*%t(G))
      t.g = as.vector(matrix(0,ng,1))
      p.val = as.vector(matrix(0,ng,1))
      se_g = as.vector(matrix(0,ng,1))
      beta_g = as.vector(beta.g)
      sigma = as.vector( sqrt((t(y)%*%t(I-S)%*%(I-S)%*%y)/ psych::tr(t(I-S)%*%(I-S))) )
      df=(u[1]^2/u[2])

      for (i in 1:ng)
      {
        se_g[i]=sigma*sqrt(gkk[i])
        t.g[i]=beta.g[i]/(sigma*sqrt(gkk[i]))
      }

      for (i in 1:ng)
      {
        if(t.g[i]>0){
          p.val[i]=2*pt(t.g[i], df=df, lower.tail = FALSE) }
        else {p.val[i]=2*pt((-1)*t.g[i], df=df, lower.tail = FALSE)}
      }


      ttabel=as.vector(matrix(qt((1-(alp/2)), df),ng,1))
      Uji.Parsial.Global=cbind(beta_g,se_g,t.g,ttabel,df,p.val=round(p.val,6))




      #membuat data frame dengan n = sebanyak amatan untuk variabel Global
      dfglob = matrix(0,length(gwr1$SDF),4*length(varglob))
      for (i in 1:length(varglob)){
        for(j in 1:length(gwr1$SDF)){
          dfglob[j,(1+(4*(i-1))):(4*i)]=Uji.Parsial.Global[i,c(1,2,3,6)]
        }
      }

      #memberi nama dfglob
      varia = c("b_","se_","t_","pv_")
      vardup1 = rep(varia,length(varglob))
      globdup = rep(varglob,each=length(varia))
      vdu1 = data.frame(vardup1,globdup)%>% tidyr::unite(nco, c(1:2), sep = "")
      dfglob = `colnames<-`(dfglob,c(vdu1$nco))

      ### Pengujian Parsial Lokal


      sigma = as.vector( sqrt((t(y)%*%t(I-S)%*%(I-S)%*%y)/ psych::tr(t(I-S)%*%(I-S))) )
      se_l=matrix(0,nl,n)
      t.hit.l=matrix(0,nl,n)
      pvalue=matrix(0,nl,n)
      ringkasan=matrix(0,n,4*nl)
      df=(u[1]^2/u[2])

      for (i in 1:n)
      {
        M=( solve(((t(xl)%*%(W[,,i]))%*%xl)) %*% t(xl) %*% (W[,,i]) %*% (I-(xg%*%G)) )
        m=diag(M%*%t(M))
        m=as.matrix(m)
        for (j in 1:nl)
        {
          se_l[j,i]=sigma*(sqrt(m[j,]))
          t.hit.l[j,i]=t(beta.l)[j,i]/(sigma*(sqrt(m[j,])))
          if(t.hit.l[j,i]>0){
            pvalue[j,i]=2*pt(t.hit.l[j,i],df=df,lower.tail=FALSE)} else
            {pvalue[j,i]=2*pt((-1)*t.hit.l[j,i],df=df,lower.tail=FALSE)}
        }
        ringkasan[i,]=t(cbind(t(beta.l)[,i],se_l[,i],t.hit.l[,i],pvalue[,i]))
      }

      ringkasan = as.data.frame(ringkasan)
      varlok2 = c("Intercept",varlok)
      vardup2 = rep(varia,length(varlok2))
      lokdup = rep(varlok2,each=length(varia))
      vdu2 = data.frame(vardup2,lokdup)%>% tidyr::unite(nco, c(1:2), sep = "")
      ringkasan = `colnames<-`(ringkasan,c(vdu2$nco)) #### INI DIUBAH JADI ringkasan

      ## Penghitungan MGWR AIC, AICc, dll
      sigmaaic = as.vector( sqrt((t(y)%*%t(I-S)%*%(I-S)%*%y)/n) )
      AICc=(2*n*log(sigmaaic))+(n*log(2*pi))+((n*((n+psych::tr(S)))/(n-2-psych::tr(S))))
      AIC=(2*n*log(sigmaaic))+(n*log(2*pi))+n+psych::tr(S)
      resid=y-y.hat
      sigu=(t(resid))%*%resid
      ym=y-mean(y)
      rsqrt1=sigu
      rsqrt2=t(ym)%*%ym
      rsqrt=1-(rsqrt1/rsqrt2) #r-squared#
      rsqrt1=rsqrt1/(n-ng-nl)
      rsqrt2=rsqrt2/(n-1)
      rbar=1-(rsqrt1/rsqrt2) #rbar-squared#
      quasiR2m = 1 - (as.vector((t(y)%*%t(I-S)%*%(I-S)%*%y))/sum(aan$`Sum Sq`))
      adjquasiR2m = 1 - ((as.vector((t(y)%*%t(I-S)%*%(I-S)%*%y))/sum(aan$`Sum Sq`))*(sum(aan$Df)/u[1]))

      #-------

      # DATAFRAME HASIL MGWR ASLI
      #dpku = data.frame(dt$Kode_Prov,dt$Provinsi,dt$Kode_Kab,dt$Kab_Kota, ringkasan, dfglob) #gaperlu udah ada


      # alhamdulillah untuk data frame variabel lokal signifikan
      grafik = data.frame()
      for (i in 1:nl){
        for (j in 1:n) {
          if(ringkasan[j,4+(4*(i-1))] < alp) {
            grafik[j,i] = colnames(xl)[i]
          } else {
            grafik[j,i] = " "
          }
        }
      }

      vdm1 = rep("sig_",length(varlok2))
      vdmu1 = data.frame(vdm1,varlok2)%>% tidyr::unite(nco, c(1:2), sep = "")
      grafik = `colnames<-`(grafik,vdmu1$nco)

      # alhamdulillah untuk data frame variabel global signifikan
      grafik_glob = data.frame()
      for (i in 1:ng){
        for (j in 1:n) {
          if(dfglob[j,4+(4*(i-1))] < alp) {
            grafik_glob[j,i] = colnames(xg)[i]
          } else {
            grafik_glob[j,i] = " "
          }
        }
      }

      vdm2 = rep("sig_",length(varglob))
      vdmu2 = data.frame(vdm2,varglob)%>% tidyr::unite(nco, c(1:2), sep = "")
      grafik_glob = `colnames<-`(grafik_glob,vdmu2$nco)

      #dpgraf = data.frame(Kode_Kab=dt$Kode_Kab, grafik_glob,grafik[,-1]) # ini gaperlu karena nanti tdk semua data frame ada kolom kodekab
      dpgraf = data.frame(grafik_glob,grafik[,-1])

      mixall = dpgraf %>%
        tidyr::unite(col = klasmgwr, c(1:(length(varlok)+length(varglob))), sep = "/")

      #dpk_2 = left_join(dpk,dpgraf, by="Kode_Kab") # ini gaperlu karena nanti tdk semua data frame ada kolom kodekab
      #dpk_3 = left_join(dpk_2,mixall, by="Kode_Kab") # ini gaperlu karena nanti tdk semua data frame ada kolom kodekab
      klasmgwr=c()
      tablemgwr = data.frame(dt, ringkasan, dfglob, grafik, grafik_glob, klasmgwr=mixall$klasmgwr)

      VMGWR = ggplot2::ggplot(dpk)+
        geom_sf(color = "white", aes(fill = mixall$klasmgwr)) +
        theme(legend.position = "left")


      z <- list(output = capture.output({
        cat("                                                 ","\n")
        cat("                                                 ","\n")
        cat("===================== Output GWR ================","\n")
        cat("                                                 ","\n")
        print(gwr1)
        cat("===================== ========== ================","\n")
        cat("                                                 ","\n")
        cat("                                                 ","\n")
        cat("                                                 ","\n")
        cat("================= Output Mixed GWR ==============","\n")
        cat("                                                 ","\n")
        cat("Varibility Test for Global and Local Variabel:","\n")
        print(LF3)
        cat("                                                 ","\n")
        cat("Global Variabel ",length(varglob)," :",varglob,"\n")
        cat("Local Variabel ",length(varlok)," :",varlok,"\n")
        cat("                                                 ","\n")
        cat("F1(GoF MGWR), F2(Global Simultaneous), F3(Local Simultaneous) test:","\n")
        print(`rownames<-`(Uji.Serentak,c("F1","F2","F3")))
        cat("                                                 ","\n")
        cat("Global Partial Test:","\n")
        print(`rownames<-`(Uji.Parsial.Global,c(varglob)))
        cat("                                                 ","\n")
        cat("AICc      = ", AICc,"\n")
        cat("AIC       = ", AIC,"\n")
        cat("R2        = ", quasiR2m,"\n")
        cat("Adj. R2   = ", rbar[1,1],"\n")
        cat("================= ================ ==============","\n")
      }),
      gwr=gwr1, Variability.Test = LF3,
      F1.F2.F3.mgwr.Test = `rownames<-`(Uji.Serentak,c("F1","F2","F3")),
      Global.Partial.Test = `rownames<-`(Uji.Parsial.Global,c(varglob)),
      table.mgwr = tablemgwr, map.mgwr = VMGWR, Global_variable=varglob,
      Local_variable=varlok, AICc=AICc, AIC=AIC, R_square=quasiR2m,
      adj_R_square=rbar[1,1])
      class(z) <- "mgwrhwc"
      print.mgwrhwc(z)
    }
   }
}




#usethis::use_data(redsb)
#utils::data()
#usethis::use_data(redsb,compress ='xz', overwrite = T)
