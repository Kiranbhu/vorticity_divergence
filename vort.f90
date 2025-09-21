program vorticity
implicit none

integer, parameter :: lon=144, lat=73, level=2
real :: uwind(lon, lat, level)
real :: vwind(lon, lat, level)
real :: vor(lon, lat, level)

integer :: i, j, l, irec
real :: dx, dy, dlat, dlon, lat_deg, lat_rad
real, parameter :: earth_radius = 6371000.0  ! Earth's radius in meters
real, parameter :: deg2rad = 3.141592653589793 / 180.0

! Open input files
open(unit=9, file='NCEP_UWND_AUG23.grd', form='unformatted', &
     access='direct', recl=lon*lat*4, status='old')
open(unit=10, file='NCEP_VWND_AUG23.grd', form='unformatted', &
     access='direct', recl=lon*lat*4, status='old')

! Open output file
open(unit=11, file='vorticity.grd', form='unformatted', &
     access='direct', recl=lon*lat*4, status='replace')

! Read wind data
irec = 1
do l = 1, level
  read(9, rec=irec) uwind(:,:,l)
  read(10, rec=irec) vwind(:,:,l)
  irec = irec + 1
end do

! Close input files
close(9)
close(10)

! Initialize vorticity array
vor = 0.0

! Calculate grid spacing (assuming 2.5 degree resolution)
dlon = 2.5  ! degrees
dlat = 2.5  ! degrees

! Compute vorticity using finite differences
do l = 1, level
  do j = 2, lat-1
    ! Convert latitude to radians for proper cosine calculation
    lat_deg = -90.0 + (j-1)*dlat  ! Latitude in degrees
    lat_rad = lat_deg * deg2rad    ! Latitude in radians
    
    ! Calculate grid spacing in meters at this latitude
    dx = earth_radius * cos(lat_rad) * dlon * deg2rad
    dy = earth_radius * dlat * deg2rad
    
    do i = 2, lon-1
      ! Vorticity formula: ζ = ∂v/∂x - ∂u/∂y
      vor(i,j,l) = ((vwind(i+1,j,l) - vwind(i-1,j,l)) / (2.0 * dx)) - &
                   ((uwind(i,j+1,l) - uwind(i,j-1,l)) / (2.0 * dy))
    end do
  end do
end do

! Handle boundaries (set to missing values)
do l = 1, level
  do j = 1, lat
    do i = 1, lon
      if (i == 1 .or. i == lon .or. j == 1 .or. j == lat) then
        vor(i,j,l) = -9.99e8  ! Missing value flag
      end if
    end do
  end do
end do

! Write vorticity data to output file
irec = 1
do l = 1, level
  write(11, rec=irec) vor(:,:,l)
  irec = irec + 1
end do

! Close output file
close(11)

print *, 'Vorticity calculation completed successfully!'
print *, 'Output saved to vorticity.grd'

end program vorticity
