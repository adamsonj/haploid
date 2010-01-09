#! /bin/sh

if [[ -e "./configure" ]]; then
    /usr/bin/autoreconf
else
    ## generate config.h.in
    /usr/bin/autoheader
    ## libtoolize
    /usr/bin/libtoolize
    /usr/bin/aclocal
    /usr/bin/automake -va
    /usr/bin/autoconf -v
fi
exit 