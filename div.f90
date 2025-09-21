program divergence
integer, parameter:: lon=144, lat=73, level=2
real:: uwind(lon, lat, level)
real:: vwind(lon, lat, level)
real:: div(lon,lat,level)
integer:: i,j,irec,jrec,l
real:: dlon, dlat, pi, re, omega
real:: lat_val, dx, dy

! Constants
pi = 3.141592653589793
re = 6371000.0    ! Earth radius in meters
omega = 7.2921150e-5  ! Earth rotation rate
dlon = 2.5        ! Longitude spacing in degrees
dlat = 2.5        ! Latitude spacing in degrees

open(9,file='NCEP_UWND_AUG23.grd', form='unformatted', access='direct', recl=lon*lat*4, status='old')
open(10,file='NCEP_VWND_AUG23.grd', form='unformatted', access='direct', recl=lon*lat*4, status='old')
open(11,file='/home/kiran/Vort_div/divergence.grd', form='unformatted', access='direct', recl=lon*lat*level*4, status='unknown')

irec = 1
do l=1, level
   read(9,rec=irec) ((uwind(i,j,l),i=1,lon),j=1,lat)
   read(10,rec=irec) ((vwind(i,j,l),i=1,lon),j=1,lat)
   irec = irec + 1
end do

! Initialize divergence with undefined value (-9.99E+08 to match your .ctl file)
do l=1,level
  do j=1, lat
    do i=1, lon
      div(i,j,l) = -9.99E+08  ! Use the same UNDEF value as in your .ctl
    end do
  end do
end do

! Compute divergence values
do l=1, level
  do j=2, lat-1
    do i=2, lon-1
      ! Calculate actual latitude in radians
      lat_val = (-90.0 + (j-1)*dlat) * pi/180.0
      
      ! Calculate grid spacing in meters
      dx = re * cos(lat_val) * (dlon * pi/180.0)
      dy = re * (dlat * pi/180.0)
      
      ! Finite difference for divergence: d(u)/dx + d(v)/dy
      if (uwind(i+1,j,l) > -9.0E+08 .and. uwind(i-1,j,l) > -9.0E+08 .and. &
          vwind(i,j+1,l) > -9.0E+08 .and. vwind(i,j-1,l) > -9.0E+08) then
        div(i,j,l) = ((uwind(i+1,j,l) - uwind(i-1,j,l)) / (2.0*dx)) + &
                     ((vwind(i,j+1,l) - vwind(i,j-1,l)) / (2.0*dy))
      end if
    end do
  end do
end do

! Write output
jrec = 1
do l=1,level
  write(11,rec=jrec) ((div(i,j,l),i=1,lon),j=1,lat)
  jrec = jrec + 1
end do

close(9)
close(10)
close(11)

print *, 'Divergence calculation completed successfully!'

end program divergence
