use super::*;

#[test]
fn names_round_trip_and_the_tools_are_named() {
    for program in Program::ALL {
        assert_eq!(Program::parse(program.name()), Ok(program));
        assert_eq!(program.to_string(), program.name());
        assert_eq!(program.name().parse::<Program>(), Ok(program));
    }
    for tool in ["makeblastdb", "makedb", "blastdbcmd"] {
        let error = Program::parse(tool).unwrap_err();
        assert!(
            error.to_string().contains("not a search program"),
            "{error}"
        );
    }
    // rpsblast and rpstblastn run on the command line only (line PS, PS11a, PS12a).
    for program in ["rpsblast", "rpstblastn"] {
        let error = Program::parse(program).unwrap_err();
        assert!(
            error.to_string().contains(&format!(
                "{program} is not available through the library API"
            )),
            "{error}"
        );
    }
    assert_eq!(
        Program::parse("nosuch"),
        Err(Error::Argv("unknown program 'nosuch'".into()))
    );
}
