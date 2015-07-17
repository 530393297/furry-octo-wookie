CC := g++

# Architecture
ARCH   := $(shell uname -m)
ifneq (,$(filter $(ARCH),x86_64 aarch64 ppc64le))
    TARGET_SIZE := 64
else ifeq ($(ARCH),armv7l)
    TARGET_SIZE := 32
else
	$(error ERROR - unsupported value $(ARCH) for ARCH!)
endif

# Operating system
OS 	:= $(shell uname -s 2>/dev/null | tr "[:upper:]" "[:lower:]")
ifeq (,$(filter $(OS),linux darwin qnx android))
    $(error ERROR - unsupported value $(OS) for OS!)
endif

# Compiler
ifeq ($(TARGET_OS),darwin)
    ifeq ($(shell expr `xcodebuild -version | grep -i xcode | awk '{print $$2}' | cut -d'.' -f1` \>= 5),1)
        HOST_COMPILER ?= clang++
    endif
endif

# Internal flags
CCFLAGS     := -std=c++11
LDFLAGS     :=

# Build flags
ifeq ($(TARGET_OS),darwin)
    CCFLAGS += -arch $(ARCH)
endif

# Common includes and paths
INCLUDES  := 
LIBRARIES :=

ifeq ($(dbg),1)
    CCFLAGS += -g 
	BUILD_TYPE := debug
else
    BUILD_TYPE := release
endif

# User flags
EXTRA_CCFLAGS := -Wall -Wextra
EXTRA_LDFLAGS :=

ALL_CCFLAGS :=
ALL_CCFLAGS += $(CCFLAGS)
ALL_CCFLAGS += $(EXTRA_CCFLAGS)

ALL_LDFLAGS :=
ALL_LDFLAGS += $(ALL_CCFLAGS)
ALL_LDFLAGS += $(LDFLAGS)
ALL_LDFLAGS += $(EXTRA_LDFLAGS)

################################################################################

# Target rules
all: build

build: lll

lll.o:lll.cpp
	$(EXEC) $(CC) $(INCLUDES) $(ALL_CCFLAGS) -o $@ -c $<

polynomial.o:polynomial.cpp
	$(EXEC) $(CC) $(INCLUDES) $(ALL_CCFLAGS) -o $@ -c $<

lll:lll.o polynomial.o
	$(EXEC) $(CC) $(ALL_LDFLAGS) -o $@ $+ $(LIBRARIES)
	$(EXEC) mkdir -p ../bin/$(ARCH)/$(OS)/$(BUILD_TYPE)
	$(EXEC) cp $@ ../bin/$(ARCH)/$(OS)/$(BUILD_TYPE)

run: build
	$(EXEC) ./lll

clean:
	rm -f lll lll.o 
	rm -rf ../bin/$(ARCH)/$(OS)/$(BUILD_TYPE)/lll

clobber: clean