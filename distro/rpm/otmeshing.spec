# norootforbuild
%global python_sitearch %{_libdir}/python%(python3 -c "import sysconfig; print(sysconfig.get_python_version())")/site-packages

Name:           otmeshing
Version:        0.0
Release:        0%{?dist}
Summary:        OpenTURNS module
Group:          System Environment/Libraries
License:        LGPLv3+
URL:            http://www.openturns.org/
Source0:        http://downloads.sourceforge.net/openturns-modules/otmeshing/otmeshing-%{version}.tar.bz2
BuildRequires:  gcc-c++, cmake, swig
BuildRequires:  openturns-devel
BuildRequires:  python3-openturns
BuildRequires:  python3-devel
BuildRequires:  CGAL-devel
BuildRequires:  cddlib-devel
BuildRequires:  qhull-devel
BuildRequires:  eigen3-devel
BuildRequires:  nanoflann-devel
BuildRequires:  tbb-devel
Requires:       libotmeshing0

%description
Template module for OpenTURNS library.

%package -n libotmeshing0
Summary:        otmeshing development files
Group:          Development/Libraries/C and C++

%description -n libotmeshing0
Dynamic libraries for otmeshing.

%package devel
Summary:        otmeshing development files
Group:          Development/Libraries/C and C++
Requires:       libotmeshing0 = %{version}
Requires:       openturns-devel

%description devel
Development files for otmeshing library.

%package -n python3-%{name}
Summary:        otmeshing library
Group:          Productivity/Scientific/Math
Requires:       python3-openturns
%description -n python3-%{name}
Python textual interface to otmeshing uncertainty library

%prep
%setup -q

%build
%cmake -DINSTALL_DESTDIR:PATH=%{buildroot} \
       -DCMAKE_SKIP_INSTALL_RPATH:BOOL=ON \
       -DCMAKE_UNITY_BUILD=ON .
%cmake_build

%install
%cmake_install

%check
export LD_LIBRARY_PATH=%{buildroot}%{_libdir}
%ctest --tests-regex pyinstallcheck --schedule-random

%post -n libotmeshing0 -p /sbin/ldconfig 
%postun -n libotmeshing0 -p /sbin/ldconfig 

%files -n libotmeshing0
%defattr(-,root,root,-)
%{_libdir}/*.so.*

%files devel
%defattr(-,root,root,-)
%dir %{_includedir}/%{name}
%{_includedir}/%{name}/*.h*
%{_includedir}/%{name}/swig/
%{_libdir}/*.so
%{_libdir}/cmake/

%files -n python3-%{name}
%defattr(-,root,root,-)
%{python_sitearch}/%{name}/
%{python_sitearch}/%{name}-*.dist-info/


%changelog
* Wed Nov 28 2012 Julien Schueller <schueller at phimeca dot com> 0.0-1
- Initial package creation

