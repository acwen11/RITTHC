CPPFLAGS += -I.
LDFLAGS += -lm -lstdc++

HDRS = finite_difference.h finite_difference.hh

SRCS = finite_difference.cc

OBJS = $(SRCS:.cc=.o)

LOCAL_SRCS = $(HDRS) $(SRCS)

all: test_finite_difference

clean:
	-rm $(LOCAL_SRCS)
	-rm test_finite_difference
	-rm finite_difference.o
	-rm test_finite_difference.dat

test_finite_difference: test_finite_difference.c $(HDRS) $(OBJS)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -std=c99 $(OBJS) test_finite_difference.c -o test_finite_difference

$(LOCAL_SRCS): %: ../src/%
	ln -s $<
