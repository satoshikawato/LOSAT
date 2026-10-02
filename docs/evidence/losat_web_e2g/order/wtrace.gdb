set pagination off
set print thread-events off
catch syscall write
commands
silent
printf "W fd=%d n=%d\n", $rdi, $rdx
continue
end
run
quit
