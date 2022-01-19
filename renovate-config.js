module.exports = {
		  branchPrefix: 'renovate/',
		  dryRun: false,
		  username: 'renovate-release',
		  gitAuthor: 'Renovate Bot <bot@renovateapp.com>',
		  onboarding: false,
		  platform: 'github',
		  requireConfig: false,
	          repositories: ["pbelmann/meta-omics-toolkit"],
		  baseBranches: ["dev"],
                  packageRules: [
                         {  "matchDatasources": ["docker"], versioning: "loose" }
		  ],
                  regexManagers: [
				{
				   "fileMatch": ["^nextflow.config$"],
				   "matchStrings": [
					"_image* = \"(?<depName>.*?):(?<currentValue>.*?)\""
				],
				   "datasourceTemplate": "docker"
				}
		]
};
