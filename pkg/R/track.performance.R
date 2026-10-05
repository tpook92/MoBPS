'#
  Authors
Torsten Pook, torsten.pook@wur.nl

Copyright (C) 2017 -- 2025  Torsten Pook

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
'#

#' Create a phenotypic transformation
#'
#' Function to perform create a transformation of phenotypes
#' @param file directory to write tracking file to
#' @param interval Interval between trackings
#' @export


track.performance <- function(file = "mobps_performance_log.tsv", interval = 5) {

  if (!grepl("^(?:[A-Za-z]:|/|\\\\\\\\)", file)) {
    file <- file.path(getwd(), file)
  }

  os <- Sys.info()[["sysname"]]

  if (os == "Windows") {

    if (requireNamespace("ps", quietly = TRUE)){

    } else{
      stop("Package 'ps' is required on Windows.")

    }



  } else {

    if(FALSE){
      pid <- Sys.getpid()

      writeLines("Timestamp\tVmHWM_kB\tVmRSS_kB\tCPU_percent_window\tCPU_percent_overall", file)

      clk_tck <- as.integer(system("getconf CLK_TCK", intern = TRUE))

      cmd <- sprintf(
        "start_ticks=$(awk '{print $14+$15}' /proc/%d/stat);
   start_time=$(date +%%s);
   prev_ticks=$start_ticks;
   prev_time=$start_time;
   while kill -0 %d 2>/dev/null; do
     sleep %d;
     now_ticks=$(awk '{print $14+$15}' /proc/%d/stat);
     now_time=$(date +%%s);
     dt=$((now_time-prev_time));
     total_dt=$((now_time-start_time));
     cpu_window=$(awk -v a=\"$now_ticks\" -v b=\"$prev_ticks\" -v dt=\"$dt\" 'BEGIN {if (dt > 0) printf \"%%.1f\",100*(a-b)/(dt*%d); else print \"NA\"}');
     cpu_overall=$(awk -v a=\"$now_ticks\" -v b=\"$start_ticks\" -v dt=\"$total_dt\" 'BEGIN {if (dt > 0) printf \"%%.1f\",100*(a-b)/(dt*%d); else print \"NA\"}');
     printf '%%s\t' \"$(date +%%F\\ %%T)\";
     awk '/^VmRSS:|^VmHWM:/ {printf \"%%s\\t\", $2}' /proc/%d/status;
     printf '%%s\t%%s\\n' \"$cpu_window\" \"$cpu_overall\";
     prev_ticks=$now_ticks;
     prev_time=$now_time;
   done >> %s & echo $!",
        pid, pid, interval, pid, clk_tck, clk_tck, pid, shQuote(file)
      )

      monitor_pid <- as.integer(system(cmd, intern = TRUE))
      return(monitor_pid)
    } else if(FALSE){

      pid <- Sys.getpid()

      writeLines(
        "Timestamp\tR_VmHWM_kB\tR_VmRSS_kB\tTree_RSS_kB\tTree_RSS_peak_kB\tR_CPU_window\tR_CPU_overall\tTree_CPU_window\tTree_CPU_overall",
        file
      )

      clk_tck <- as.integer(system("getconf CLK_TCK", intern = TRUE))

      # Build the shell script as a vector of lines (real quotes, no R-escaping needed)
      script_lines <- c(
        'root=%d',
        'tree_peak_rss=0',
        'get_pids() {',
        '  ps -eo pid=,ppid= | awk -v root=$root \'',
        '    BEGIN { pids[root]=1; changed=1 }',
        '    { all_pid[$1]=$1; parent[$1]=$2 }',
        '    END {',
        '      while (changed) {',
        '        changed=0',
        '        for (pid in all_pid) {',
        '          if (pids[parent[pid]] && !pids[pid]) {',
        '            pids[pid]=1',
        '            changed=1',
        '          }',
        '        }',
        '      }',
        '      for (pid in pids) print pid',
        '    }\'',
        '}',
        'sum_ticks() {',
        '  total=0',
        '  for p in $(get_pids); do',
        '    if [ -r /proc/$p/stat ]; then',
        '      t=$(awk \'{print $14+$15}\' /proc/$p/stat 2>/dev/null)',
        '      total=$((total + ${t:-0}))',
        '    fi',
        '  done',
        '  echo $total',
        '}',
        'sum_rss() {',
        '  total=0',
        '  for p in $(get_pids); do',
        '    if [ -r /proc/$p/status ]; then',
        '      r=$(awk \'/^VmRSS:/ {print $2}\' /proc/$p/status 2>/dev/null)',
        '      total=$((total + ${r:-0}))',
        '    fi',
        '  done',
        '  echo $total',
        '}',
        'r_start_ticks=$(awk \'{print $14+$15}\' /proc/$root/stat)',
        'tree_start_ticks=$(sum_ticks)',
        'start_time=$(date +%%s)',
        'r_prev_ticks=$r_start_ticks',
        'tree_prev_ticks=$tree_start_ticks',
        'prev_time=$start_time',
        'while kill -0 $root 2>/dev/null; do',
        '  sleep %d',
        '  now_time=$(date +%%s)',
        '  dt=$((now_time-prev_time))',
        '  total_dt=$((now_time-start_time))',
        '  r_now_ticks=$(awk \'{print $14+$15}\' /proc/$root/stat)',
        '  tree_now_ticks=$(sum_ticks)',
        '  r_cpu_window=$(awk -v a="$r_now_ticks" -v b="$r_prev_ticks" -v dt="$dt" -v clk=%d \'BEGIN {if (dt > 0) printf "%%.1f", 100*(a-b)/(dt*clk); else print "NA"}\')',
        '  r_cpu_overall=$(awk -v a="$r_now_ticks" -v b="$r_start_ticks" -v dt="$total_dt" -v clk=%d \'BEGIN {if (dt > 0) printf "%%.1f", 100*(a-b)/(dt*clk); else print "NA"}\')',
        '  tree_cpu_window=$(awk -v a="$tree_now_ticks" -v b="$tree_prev_ticks" -v dt="$dt" -v clk=%d \'BEGIN {if (dt > 0) printf "%%.1f", 100*(a-b)/(dt*clk); else print "NA"}\')',
        '  tree_cpu_overall=$(awk -v a="$tree_now_ticks" -v b="$tree_start_ticks" -v dt="$total_dt" -v clk=%d \'BEGIN {if (dt > 0) printf "%%.1f", 100*(a-b)/(dt*clk); else print "NA"}\')',
        '  r_hwm=$(awk \'/^VmHWM:/ {print $2}\' /proc/$root/status)',
        '  r_rss=$(awk \'/^VmRSS:/ {print $2}\' /proc/$root/status)',
        '  tree_rss=$(sum_rss)',
        '  if [ $tree_rss -gt $tree_peak_rss ]; then',
        '    tree_peak_rss=$tree_rss',
        '  fi',
        '  printf "%%s\\t%%s\\t%%s\\t%%s\\t%%s\\t%%s\\t%%s\\t%%s\\t%%s\\n" \\',
        '    "$(date +%%F\\ %%T)" \\',
        '    "$r_hwm" \\',
        '    "$r_rss" \\',
        '    "$tree_rss" \\',
        '    "$tree_peak_rss" \\',
        '    "$r_cpu_window" \\',
        '    "$r_cpu_overall" \\',
        '    "$tree_cpu_window" \\',
        '    "$tree_cpu_overall"',
        '  r_prev_ticks=$r_now_ticks',
        '  tree_prev_ticks=$tree_now_ticks',
        '  prev_time=$now_time',
        'done >> %s &',
        'echo $!'
      )

      # Fill in the placeholders (%d/%s). %% -> literal % survives here.
      script_lines <- sprintf(
        paste(script_lines, collapse = "\n"),
        pid,        # root=%d
        interval,   # sleep %d
        clk_tck, clk_tck, clk_tck, clk_tck,  # the four -v clk=%d
        shQuote(file)  # done >> %s &
      )

      script_file <- tempfile(fileext = ".sh")
      writeLines(script_lines, script_file)
      Sys.chmod(script_file, "0700")

      monitor_pid <- as.integer(system(
        paste("bash", shQuote(script_file)),
        intern = TRUE
      ))
      return(monitor_pid)

    } else{

      root_pid <- Sys.getpid()

      # Capture start time of the R process itself, to detect PID reuse later
      root_starttime <- system(
        sprintf("awk '{print $22}' /proc/%d/stat", root_pid),
        intern = TRUE
      )

      # Cgroup path (whole session: R + all system() children)
      cg_path <- system("awk -F: '{print $3}' /proc/self/cgroup", intern = TRUE)[1]
      cg_dir  <- paste0("/sys/fs/cgroup", cg_path)

      cpu_stat_file <- file.path(cg_dir, "cpu.stat")
      mem_cur_file  <- file.path(cg_dir, "memory.current")

      if (!dir.exists(cg_dir)) {
        stop("cgroup directory not found: ", cg_dir)
      }

      writeLines(
        "Timestamp\tSession_current_kB\tSession_peak_kB\tR_current_kB\tR_peak_kB\tCPU_usage_pct_window\tCPU_usage_pct_overall",
        file
      )

      script_lines <- c(
        "root=%d",
        "root_starttime=%s",
        "cpu_stat_file=%s",
        "mem_cur_file=%s",

        "get_usage_usec() { awk '/^usage_usec/ {print $2}' \"$cpu_stat_file\"; }",
        "get_mem_cur()    { cat \"$mem_cur_file\"; }",

        "root_alive() {",
        "  [ -r /proc/$root/stat ] || return 1",
        "  cur=$(awk '{print $22}' /proc/$root/stat 2>/dev/null)",
        "  [ \"$cur\" = \"$root_starttime\" ]",
        "}",

        "(",
        "session_peak=0",
        "start_usec=$(get_usage_usec)",
        "prev_usec=$start_usec",
        "start_time=$(date +%%s.%%N)",
        "prev_time=$start_time",

        "while root_alive; do",
        "  sleep %f",
        "  now_time=$(date +%%s.%%N)",
        "  now_usec=$(get_usage_usec)",

        "  session_cur=$(get_mem_cur)",
        "  session_cur_kb=$((session_cur / 1024))",
        "  if [ \"$session_cur_kb\" -gt \"$session_peak\" ]; then session_peak=$session_cur_kb; fi",

        "  r_status=/proc/$root/status",
        "  if [ -r \"$r_status\" ]; then",
        "    r_cur_kb=$(awk '/^VmRSS:/ {print $2}' \"$r_status\")",
        "    r_peak_kb=$(awk '/^VmHWM:/ {print $2}' \"$r_status\")",
        "  else",
        "    r_cur_kb=NA",
        "    r_peak_kb=NA",
        "  fi",

        "  dt=$(awk -v a=\"$now_time\" -v b=\"$prev_time\" 'BEGIN {print a-b}')",
        "  total_dt=$(awk -v a=\"$now_time\" -v b=\"$start_time\" 'BEGIN {print a-b}')",

        "  cpu_pct_window=$(awk -v a=\"$now_usec\" -v b=\"$prev_usec\" -v dt=\"$dt\" 'BEGIN {if (dt > 0) printf \"%%.1f\", 100*(a-b)/(dt*1000000); else print \"NA\"}')",
        "  cpu_pct_overall=$(awk -v a=\"$now_usec\" -v b=\"$start_usec\" -v dt=\"$total_dt\" 'BEGIN {if (dt > 0) printf \"%%.1f\", 100*(a-b)/(dt*1000000); else print \"NA\"}')",

        "  printf '%%s\\t%%s\\t%%s\\t%%s\\t%%s\\t%%s\\t%%s\\n' \"$(date +'%%F %%T')\" \"$session_cur_kb\" \"$session_peak\" \"$r_cur_kb\" \"$r_peak_kb\" \"$cpu_pct_window\" \"$cpu_pct_overall\" >> %s",

        "  prev_usec=$now_usec",
        "  prev_time=$now_time",
        "done",
        ") < /dev/null > /dev/null 2>&1 &",
        "disown",
        "echo $!"
      )

      script_txt <- sprintf(
        paste(script_lines, collapse = "\n"),
        root_pid,                 # root=%d
        root_starttime,           # root_starttime=%s
        shQuote(cpu_stat_file),   # cpu_stat_file=%s
        shQuote(mem_cur_file),    # mem_cur_file=%s
        interval,                 # sleep %f
        shQuote(file)             # >> %s
      )

      script_file <- tempfile(fileext = ".sh")
      writeLines(script_txt, script_file)
      Sys.chmod(script_file, "0700")

      monitor_pid <- as.integer(system(
        paste("bash", shQuote(script_file)),
        intern = TRUE
      ))


      return(monitor_pid)

      to_kill = sprintf("kill %d", monitor_pid)
      system(to_kill)
    }

  }
}
