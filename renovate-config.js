module.exports = {
	  branchPrefix: 'renovate/',
	  dryRun: true,
	  username: 'renovate-release',
	  gitAuthor: 'Renovate Bot',
	  onboarding: true,
	  platform: 'github',
	  packageRules: [{
		"regexManagers": [
		   {
		      "fileMatch": ["^nextflow.config$"],
		      "matchStrings": [
		              "_image* = \"(?<depName>.*?):(?<currentValue>.*?)\""
		            ],
		      "datasourceTemplate": "docker"
		    }
		    ]
		  }
	  ],};
