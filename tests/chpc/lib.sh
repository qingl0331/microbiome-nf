# Assertion helpers shared by the CHPC test scripts
FAILS=0
pass()    { echo "  PASS  $*"; }
fail()    { echo "  FAIL  $*"; FAILS=$((FAILS + 1)); }
exists()  { [ -e "$1" ] && pass "exists $1" || fail "missing $1"; }
nonempty(){ [ -s "$1" ] && pass "non-empty ${1#$OUT/}" || fail "missing/empty ${1#$OUT/}"; }
absent()  { [ ! -e "$1" ] && pass "absent ${1#$OUT/}" || fail "should not exist: ${1#$OUT/}"; }
qza_ok()  { unzip -tq "$1" >/dev/null 2>&1 && pass "valid archive $(basename "$1")" || fail "corrupt/empty archive $1"; }
summary() {
    echo
    if [ "$FAILS" -eq 0 ]; then echo "ALL CHECKS PASSED"; else echo "$FAILS CHECK(S) FAILED"; fi
    exit "$FAILS"
}
# expect_fail "<description>" "<grep pattern for the error message>" nextflow run ...
expect_fail() {
    local desc=$1 pat=$2; shift 2
    if "$@" >"$TEST_ROOT/expect_fail.log" 2>&1; then
        fail "$desc: run succeeded but should have failed"
    elif [ -n "$pat" ] && ! grep -q -- "$pat" "$TEST_ROOT/expect_fail.log"; then
        fail "$desc: failed, but without '$pat' (see $TEST_ROOT/expect_fail.log)"
    else
        pass "$desc"
    fi
}
# status of every task in the last run: CACHED / COMPLETED / FAILED
task_status() { nextflow log last -f name,status 2>/dev/null; }
