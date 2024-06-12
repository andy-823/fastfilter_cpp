

# Uncomment one of the following to switch between debug and opt mode
# OPT = -O3 -DNDEBUG -std=c++17
OPT = -O3 -std=c++17 -g -ggdb -fno-omit-frame-pointer -Wextra -fpermissive # -fsanitize=undefined -fsanitize=address

CXXFLAGS += -fno-strict-aliasing -Wall -std=c++11 -I. -I../src/ \
    -I../src/bloom/ -I../src/cuckoo/ -I../src/gcs \
    -I../src/gqf/ -I../src/morton/ -I../src/xorfilter -I../src/ribbon \
    -I../src/vqf \
    $(OPT)

UNAME_P := $(shell uname -p)
ifeq ($(UNAME_P),x86_64)
        CXXFLAGS += -march=native
else
        CXXFLAGS +=
endif
LDFLAGS = -Wall -Wextra

HEADERS = $(wildcard ../src/*.h \
    ../src/bloom/*.h ../src/cuckoo/*.h ../src/gcs/*.h \
    ../src/gqf/*.h ../src/morton/*.h ../src/xorfilter/*.h ../src/ribbon/*.h \
    ../src/vqf/*.h \
    ) *.h

.PHONY: all

BINS = bulk-insert-and-query.exe stream.exe construction-failure.exe \
       test-block-count.exe test-block-cout2.exe test-block-count3.exe \
       compare-sizes.exe

all: $(BINS)

clean:
	/bin/rm -f $(BINS)

%.exe: %.cc ${HEADERS}  Makefile
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)
